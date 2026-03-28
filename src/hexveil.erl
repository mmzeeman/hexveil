%% @doc HexVeil - A privacy-preserving hierarchical location encoding
%% system using a flat hexagonal grid with axial (q, r) coordinates.
%%
%% Optimized for "Human Scale" landmarks (~2.4m) and "Smooth Privacy" steps.
%% This implementation uses 24 levels of precision.
%%
%% Cell sizes (flat-to-flat, pointy-top hex):
%%   Level 24: ~2.4m       — Meeting point / Bench (finest)
%%   Level 22: ~9.6m       — Building entrance
%%   Level 20: ~38.4m      — Large building / Shop
%%   Level 18: ~153.6m     — Privacy 1: Building block
%%   Level 17: ~307.2m     — Privacy 2: Neighborhood
%%   Level 16: ~614.4m     — Privacy 3: District
%%   Level 14: ~2.46km     — Local region
%%   Level 12: ~9.8km      — City Scale
%%   Level 8:  ~157km      — Metropolitan area
%%   Level 4:  ~2,516km     — Sub-continental
%%   Level 1:  ~20,132km    — Global Scale
%%
%% Each level is exactly 2 times coarser (linearly) than the one below it.
%% This is an Aperture 4 hierarchy (Area x 4 per level). No rotation.
%%
%% Code format:
%%   {Q, R} integers at level 24 (finest).
%%   display/2 renders as an interleaved Hexadecimal string (0-F).
%%   Each character represents TWO levels (4 bits: Qn Rn Qn+1 Rn+1).
%%
%% Usage:
%%   Code  = hexveil:encode(52.3616, 4.8784).
%%   S     = hexveil:display(Code).              %% "6CEC50ECC574"
%%
%% Examples (Level 24):
%%   Vondelpark Entrance : 6CEC50ECC574
%%   Leidseplein, Ams.   : 6CEC50EDAE5F
%%   De Dam, Amsterdam   : 6CEC50FBAC1E
%%   De Dom, Utrecht     : 6CEC632CA909
%%
%% Hierarchical Examples:
%%   Location             | L22 (~9.6m) | L18 (~154m) | L17 (~307m) | L16 (~614m) | L12 (~9.8km)
%%   ---------------------|-------------|-------------|-------------|-------------|------------
%%   Leidseplein          | 6CEC50EDAE5 | 6CEC50EDA   | 6CEC50ED2   | 6CEC50ED    | 6CEC50
%%   Vondelpark Entrance  | 6CEC50ECC57 | 6CEC50ECC   | 6CEC50EC3   | 6CEC50EC    | 6CEC50
%%   De Dam, Amsterdam    | 6CEC50FBAC1 | 6CEC50FBA   | 6CEC50FB2   | 6CEC50FB    | 6CEC50
%%   De Dom, Utrecht      | 6CEC632CA90 | 6CEC632CA   | 6CEC632C2   | 6CEC632C    | 6CEC63
%%
%% Projection: pointy-top hex, equirectangular, centred on (0.0, 0.0).
%% Coverage: Global.

-module(hexveil).

-export([
    encode/2,
    decode/1,
    are_nearby/3,
    coarsen/2,
    neighbors/1,
    cell_bounds/1,
    cell_geometry/1,
    cell_geometry/2,
    display/1,
    display/2,
    parse/1
]).

-type code()  :: {non_neg_integer(), non_neg_integer()}.
-type level() :: 1..24.
-export_type([code/0, level/0]).

%% ---------------------------------------------------------------------------
%% Constants
%% ---------------------------------------------------------------------------

-define(MAX_LEVEL, 24).
-define(R, 1.385641).
-define(BQ_X, 2.4).
-define(BR_X, 1.2).
-define(BR_Y, 2.078461).
-define(REF_LAT,  0.0).
-define(REF_LON,  0.0).
-define(M_PER_DEG_LAT, 111319.49).
-define(M_PER_DEG_LON, 111319.49).

%% Offset of 2^23 centers the 24-bit grid globally (~40,265km range).
-define(Q_OFF, 8388608).
-define(R_OFF, 8388608).
-define(DIRECTIONS, [{1,0},{0,1},{-1,1},{-1,0},{0,-1},{1,-1}]).

%% ---------------------------------------------------------------------------
%% Public API
%% ---------------------------------------------------------------------------

encode(Lat, Lon) ->
    {X, Y}   = latlon_to_xy(Lat, Lon),
    {Qf, Rf} = xy_to_axial(X, Y),
    {Q,  R}  = hex_round(Qf, Rf),
    {Q + ?Q_OFF, R + ?R_OFF}.

decode({Q, R}) ->
    {X, Y} = axial_to_xy(Q - ?Q_OFF + 0.5, R - ?R_OFF + 0.5),
    xy_to_latlon(X, Y).

are_nearby({Q1, R1}, {Q2, R2}, Level) ->
    Shift = ?MAX_LEVEL - Level,
    Q1 bsr Shift =:= Q2 bsr Shift andalso
    R1 bsr Shift =:= R2 bsr Shift.

coarsen({Q, R}, Level) ->
    Shift = ?MAX_LEVEL - Level,
    {Q bsr Shift, R bsr Shift}.

neighbors({Q, R}) ->
    [{Q + DQ, R + DR} || {DQ, DR} <- ?DIRECTIONS].

cell_bounds(Code) ->
    {CLat, CLon} = decode(Code),
    Half    = 1.2,
    HalfLat = Half / ?M_PER_DEG_LAT,
    HalfLon = Half / ?M_PER_DEG_LON,
    {CLat - HalfLat, CLon - HalfLon,
     CLat + HalfLat, CLon + HalfLon}.

%% @doc Return the 6 corner coordinates {Lat, Lon} of a cell at level 24.
-spec cell_geometry(code()) -> [{float(), float()}].
cell_geometry(Code) -> cell_geometry(Code, ?MAX_LEVEL).

%% @doc Return the 6 corner coordinates {Lat, Lon} of a cell at the given level.
-spec cell_geometry(code(), level()) -> [{float(), float()}].
cell_geometry({Q, R}, Level) ->
    {CX, CY} = axial_to_xy(Q - ?Q_OFF + 0.5, R - ?R_OFF + 0.5),
    BaseRadius = ?BQ_X / math:sqrt(3.0),
    Radius = BaseRadius * math:pow(2, ?MAX_LEVEL - Level),
    Angles = [30, 90, 150, 210, 270, 330],
    [xy_to_latlon(
        CX + Radius * math:cos(A * math:pi() / 180),
        CY + Radius * math:sin(A * math:pi() / 180)
     ) || A <- Angles].

display(Code) -> display(Code, ?MAX_LEVEL).

display({Q, R}, Level) ->
    %% Interleave the 24 bits of the encoded (offset) coordinates.
    Bits = interleave_bits(<<Q:24>>, <<R:24>>, <<>>),
    Required = Level * 2,
    <<Prefix:Required/bits, _/bits>> = Bits,
    bits_to_hex(Prefix, <<>>).

parse(S) ->
    Bits = hex_to_bits(S, <<>>),
    Level = bit_size(Bits) div 2,
    {Qp, Rp} = deinterleave_bits(Bits, 0, 0),
    %% Return the level-24 representation of the coarsened cell.
    Shift = ?MAX_LEVEL - Level,
    {Qp bsl Shift, Rp bsl Shift}.

%% ---------------------------------------------------------------------------
%% Bit Interleaving (internal)
%% ---------------------------------------------------------------------------

interleave_bits(<<Q:1, QRest/bits>>, <<R:1, RRest/bits>>, Acc) ->
    interleave_bits(QRest, RRest, <<Acc/bits, Q:1, R:1>>);
interleave_bits(<<>>, <<>>, Acc) -> Acc.

bits_to_hex(<<V:4, Rest/bits>>, Acc) ->
    bits_to_hex(Rest, <<Acc/binary, (hex_char(V))/binary>>);
bits_to_hex(<<V:2, Rest/bits>>, Acc) ->
    bits_to_hex(Rest, <<Acc/binary, (hex_char(V))/binary>>);
bits_to_hex(<<>>, Acc) -> Acc.

hex_char(V) when V < 10 -> <<($0 + V)>>;
hex_char(V) -> <<($A + V - 10)>>.

hex_to_bits(<<C:1/binary, Rest/binary>>, Acc) when byte_size(Rest) > 0 ->
    Val = hex_val(C),
    hex_to_bits(Rest, <<Acc/bits, Val:4>>);
hex_to_bits(<<C:1/binary>>, Acc) ->
    Val = hex_val(C),
    %% In our system, odd levels result in 2-bit residual characters.
    if Val > 3 -> <<Acc/bits, Val:4>>;
       true    -> <<Acc/bits, Val:2>>
    end;
hex_to_bits(<<>>, Acc) -> Acc.

hex_val(<<C>>) when C >= $0, C =< $9 -> C - $0;
hex_val(<<C>>) when C >= $A, C =< $F -> C - $A + 10;
hex_val(<<C>>) when C >= $a, C =< $f -> C - $a + 10.

deinterleave_bits(<<QBit:1, RBit:1, Rest/bits>>, Q, R) ->
    deinterleave_bits(Rest, (Q bsl 1) bor QBit, (R bsl 1) bor RBit);
deinterleave_bits(<<>>, Q, R) -> {Q, R}.

%% ---------------------------------------------------------------------------
%% Geometry (internal)
%% ---------------------------------------------------------------------------

latlon_to_xy(Lat, Lon) ->
    {(Lon - ?REF_LON) * ?M_PER_DEG_LON,
     (Lat - ?REF_LAT) * ?M_PER_DEG_LAT}.

xy_to_latlon(X, Y) ->
    {?REF_LAT + Y / ?M_PER_DEG_LAT,
     ?REF_LON + X / ?M_PER_DEG_LON}.

xy_to_axial(X, Y) ->
    Rf = Y / (1.5 * ?R),
    Qf = (X - ?BR_X * Rf) / ?BQ_X,
    {Qf, Rf}.

axial_to_xy(Q, R) ->
    {?BQ_X * Q + ?BR_X * R,
     ?BR_Y * R}.

hex_round(Qf, Rf) ->
    Sf = -Qf - Rf,
    Rq = round(Qf), Rr = round(Rf), Rs = round(Sf),
    Dq = abs(Rq - Qf), Dr = abs(Rr - Rf), Ds = abs(Rs - Sf),
    if
        Dq > Dr andalso Dq > Ds -> {-Rr - Rs, Rr};
        Dr > Ds                  -> {Rq, -Rq - Rs};
        true                     -> {Rq, Rr}
    end.
