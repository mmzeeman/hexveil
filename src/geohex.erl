%% @doc GeoHex - A privacy-preserving hierarchical location encoding
%% system using a flat hexagonal grid with axial (q, r) coordinates.
%%
%% Optimized for "City Scale" use cases (lending, coffee, local events).
%% This implementation uses 7 levels of precision with a ~150m base cell.
%%
%% Cell sizes (flat-to-flat, pointy-top hex):
%%   Level 7: ~150m      — building/complex scale (high precision)
%%   Level 6: ~1.05km    — neighbourhood          (walking distance)
%%   Level 5: ~7.35km    — city district          (10-min drive)
%%   Level 4: ~51.45km   — city/metro             (covers major cities)
%%   Level 3: ~360.15km  — regional               (multi-city)
%%   Level 2: ~2,521km   — sub-continental        (multi-state)
%%   Level 1: ~17,647km  — global sector          (hemisphere scale)
%%
%% Each level is exactly 7 times coarser (linearly) than the one below it.
%% This is an Aperture 49 hierarchy (Area x 49 per level).
%%
%% Code format:
%%   {Q, R} integers at level 7 (finest).
%%   display/2 renders as a Base-49 string.
%%   Each character represents one level (interleaving Q and R digits).
%%   Level 7 string is 7 characters, Level 4 is 4 characters.
%%   This enables efficient database prefix queries (e.g. WHERE hex LIKE 'abcd%').
%%
%% Projection: pointy-top hex, equirectangular, centred on (0.0, 0.0).
%% Coverage: Global.
%%
%% Usage:
%%   Code  = geohex:encode(52.3026, 4.6889).
%%   true  = geohex:are_nearby(Code, Other, 6). %% same ~1km area?
%%   Nbrs  = geohex:neighbors(Code).
%%   S     = geohex:display(Code).              %% 7-char Base-49 string
%%   S4    = geohex:display(Code, 4).           %% 4-char city-level string
%%   Code2 = geohex:parse(S).

-module(geohex).

-export([
    encode/2,
    decode/1,
    are_nearby/3,
    coarsen/2,
    neighbors/1,
    cell_bounds/1,
    display/1,
    display/2,
    parse/1
]).

-type code()  :: {non_neg_integer(), non_neg_integer()}.
-type level() :: 1..7.
-export_type([code/0, level/0]).

%% ---------------------------------------------------------------------------
%% Constants
%% ---------------------------------------------------------------------------

%% Finest level is 7.
-define(MAX_LEVEL, 7).

%% Circumradius R at level 7: flat-to-flat = sqrt(3)*R = 150m
-define(R, 86.602540).

%% Axial basis vectors in metres (pointy-top hex)
-define(BQ_X, 150.0).      %% east  metres per q step  = sqrt(3)*R
-define(BR_X, 75.0).       %% east  metres per r step  = sqrt(3)/2*R
-define(BR_Y, 129.90381).  %% north metres per r step  = 3/2*R

%% Projection reference point
-define(REF_LAT,  0.0).
-define(REF_LON,  0.0).
-define(M_PER_DEG_LAT, 111319.49).
-define(M_PER_DEG_LON, 111319.49).

%% Offset so all level-7 encoded values are non-negative.
%% 7^7 = 823,543 cells per axis.
-define(Q_OFF, 411771).
-define(R_OFF, 411771).

%% Alphabet for Base-49: 0-9, a-z, A-M (49 total)
-define(ALPHABET, <<"0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLM">>).

%% The 6 axial neighbour directions for a pointy-top hex grid.
-define(DIRECTIONS, [{1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1}]).

%% ---------------------------------------------------------------------------
%% Public API
%% ---------------------------------------------------------------------------

%% @doc Encode a lat/lon pair into a {Q, R} code at level 7 (finest).
-spec encode(float(), float()) -> code().
encode(Lat, Lon) ->
    {X, Y}   = latlon_to_xy(Lat, Lon),
    {Qf, Rf} = xy_to_axial(X, Y),
    {Q,  R}  = hex_round(Qf, Rf),
    {Q + ?Q_OFF, R + ?R_OFF}.

%% @doc Decode a code back to {Lat, Lon} at the centre of its level-7 cell.
-spec decode(code()) -> {float(), float()}.
decode({Q, R}) ->
    {X, Y} = axial_to_xy(Q - ?Q_OFF + 0.5, R - ?R_OFF + 0.5),
    xy_to_latlon(X, Y).

%% @doc True if two codes refer to the same cell at the given level.
-spec are_nearby(code(), code(), level()) -> boolean().
are_nearby({Q1, R1}, {Q2, R2}, Level) ->
    Div = pow7(?MAX_LEVEL - Level),
    Q1 div Div =:= Q2 div Div andalso
    R1 div Div =:= R2 div Div.

%% @doc Coarsen a code to the given level by integer division.
-spec coarsen(code(), level()) -> code().
coarsen({Q, R}, Level) ->
    Div = pow7(?MAX_LEVEL - Level),
    {Q div Div, R div Div}.

%% @doc Return the 6 neighbouring codes at level 7.
-spec neighbors(code()) -> [code()].
neighbors({Q, R}) ->
    [{Q + DQ, R + DR} || {DQ, DR} <- ?DIRECTIONS].

%% @doc Approximate bounding box {MinLat, MinLon, MaxLat, MaxLon}.
-spec cell_bounds(code()) -> {float(), float(), float(), float()}.
cell_bounds(Code) ->
    {CLat, CLon} = decode(Code),
    Half    = 75.0,
    HalfLat = Half / ?M_PER_DEG_LAT,
    HalfLon = Half / ?M_PER_DEG_LON,
    {CLat - HalfLat, CLon - HalfLon,
     CLat + HalfLat, CLon + HalfLon}.

%% @doc Render a code as a human-readable Base-49 binary at its finest level.
%% Produces a 7-character string.
-spec display(code()) -> binary().
display(Code) -> display(Code, ?MAX_LEVEL).

%% @doc Render a code at the given level.
%% Produces a binary of Level characters.
-spec display(code(), level()) -> binary().
display({Q, R}, Level) ->
    Div = pow7(?MAX_LEVEL - Level),
    QCoarse = Q div Div,
    RCoarse = R div Div,
    encode_base49(QCoarse, RCoarse, Level).

%% @doc Parse a display binary back into a level-7 {Q, R} code.
-spec parse(binary()) -> code().
parse(S) ->
    Level = byte_size(S),
    {Qp, Rp} = decode_base49(S, Level),
    Div = pow7(?MAX_LEVEL - Level),
    {Qp * Div, Rp * Div}.

%% ---------------------------------------------------------------------------
%% Base-49 Logic (internal)
%% ---------------------------------------------------------------------------

-spec encode_base49(non_neg_integer(), non_neg_integer(), level()) -> binary().
encode_base49(Q, R, Level) ->
    QBin = pad7(integer_to_binary(Q, 7), Level),
    RBin = pad7(integer_to_binary(R, 7), Level),
    map_to_alphabet(QBin, RBin, <<>>).

map_to_alphabet(<<Q:1/binary, QRest/binary>>, <<R:1/binary, RRest/binary>>, Acc) ->
    Index = (binary_to_integer(Q, 7) * 7) + binary_to_integer(R, 7),
    Char  = binary:part(?ALPHABET, Index, 1),
    map_to_alphabet(QRest, RRest, <<Acc/binary, Char/binary>>);
map_to_alphabet(<<>>, <<>>, Acc) ->
    Acc.

-spec decode_base49(binary(), level()) -> {non_neg_integer(), non_neg_integer()}.
decode_base49(S, _Level) ->
    {QBin, RBin} = map_from_alphabet(S, <<>>, <<>>),
    {binary_to_integer(QBin, 7), binary_to_integer(RBin, 7)}.

map_from_alphabet(<<Char:1/binary, Rest/binary>>, QAcc, RAcc) ->
    {Index, 1} = binary:match(?ALPHABET, Char),
    Q = integer_to_binary(Index div 7, 7),
    R = integer_to_binary(Index rem 7, 7),
    map_from_alphabet(Rest, <<QAcc/binary, Q/binary>>, <<RAcc/binary, R/binary>>);
map_from_alphabet(<<>>, QAcc, RAcc) ->
    {QAcc, RAcc}.

%% ---------------------------------------------------------------------------
%% Helpers
%% ---------------------------------------------------------------------------

-spec pow7(0..6) -> pos_integer().
pow7(0) -> 1;
pow7(1) -> 7;
pow7(2) -> 49;
pow7(3) -> 343;
pow7(4) -> 2401;
pow7(5) -> 16807;
pow7(6) -> 117649.

-spec pad7(binary(), pos_integer()) -> binary().
pad7(S, Width) ->
    Padding = binary:copy(<<"0">>, max(0, Width - byte_size(S))),
    <<Padding/binary, S/binary>>.

-spec latlon_to_xy(float(), float()) -> {float(), float()}.
latlon_to_xy(Lat, Lon) ->
    {(Lon - ?REF_LON) * ?M_PER_DEG_LON,
     (Lat - ?REF_LAT) * ?M_PER_DEG_LAT}.

-spec xy_to_latlon(float(), float()) -> {float(), float()}.
xy_to_latlon(X, Y) ->
    {?REF_LAT + Y / ?M_PER_DEG_LAT,
     ?REF_LON + X / ?M_PER_DEG_LON}.

-spec xy_to_axial(float(), float()) -> {float(), float()}.
xy_to_axial(X, Y) ->
    Rf = Y / (1.5 * ?R),
    Qf = (X - ?BR_X * Rf) / ?BQ_X,
    {Qf, Rf}.

-spec axial_to_xy(float(), float()) -> {float(), float()}.
axial_to_xy(Q, R) ->
    {?BQ_X * Q + ?BR_X * R,
     ?BR_Y * R}.

-spec hex_round(float(), float()) -> {integer(), integer()}.
hex_round(Qf, Rf) ->
    Sf = -Qf - Rf,
    Rq = round(Qf), Rr = round(Rf), Rs = round(Sf),
    Dq = abs(Rq - Qf), Dr = abs(Rr - Rf), Ds = abs(Rs - Sf),
    if
        Dq > Dr andalso Dq > Ds -> {-Rr - Rs, Rr};
        Dr > Ds                  -> {Rq, -Rq - Rs};
        true                     -> {Rq, Rr}
    end.
