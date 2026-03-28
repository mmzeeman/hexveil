%% @doc HexVeil - A privacy-preserving hierarchical location encoding
%% system using a flat hexagonal grid with axial (q, r) coordinates.
%%
%% Optimized for "Human Scale" landmarks (~2.4m) and "Smooth Privacy" steps.
%% This implementation uses 24 levels of precision.
%%
%% Each level is exactly 2 times coarser (linearly) than the one below it.
%% This is an Aperture 4 hierarchy (Area x 4 per level). No rotation.
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
%% Projection: pointy-top hex, sinusoidal (Lon scaled by cos(Lat)).
%% Coverage: Global.

-module(hexveil).

-export([
    encode/2,
    decode/1,
    decode/2,
    are_nearby/3,
    coarsen/2,
    neighbors/1,
    cell_bounds/1,
    cell_bounds/2,
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

%% Offset centers the 24-bit grid globally (~40,265km range).
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

decode({Q, R}) -> decode({Q, R}, ?MAX_LEVEL).

%% @doc Decode a code to {Lat, Lon} at the centre of its cell at the given level.
-spec decode(code(), level()) -> {float(), float()}.
decode({Q, R}, Level) ->
    Shift = ?MAX_LEVEL - Level,
    {Qc, Rc} = if Shift == 0 -> {Q - ?Q_OFF + 0.5, R - ?R_OFF + 0.5};
                  true ->
                      Q_base = (Q bsr Shift) bsl Shift,
                      R_base = (R bsr Shift) bsl Shift,
                      Center = 1 bsl (Shift - 1),
                      {Q_base - ?Q_OFF + Center, R_base - ?R_OFF + Center}
               end,
    {X, Y} = axial_to_xy(Qc, Rc),
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

cell_bounds(Code) -> cell_bounds(Code, ?MAX_LEVEL).

%% @doc Approximate bounding box {MinLat, MinLon, MaxLat, MaxLon} at the given level.
-spec cell_bounds(code(), level()) -> {float(), float(), float(), float()}.
cell_bounds(Code, Level) ->
    {CLat, CLon} = decode(Code, Level),
    %% Radius at Level 24 is 1.2m. Scale by 2x per level.
    Half = 1.2 * math:pow(2, ?MAX_LEVEL - Level),
    HalfLat = Half / ?M_PER_DEG_LAT,
    HalfLon = Half / (?M_PER_DEG_LAT * math:cos(CLat * math:pi() / 180.0)),
    {CLat - HalfLat, CLon - HalfLon,
     CLat + HalfLat, CLon + HalfLon}.

%% @doc Return the 6 corner coordinates {Lat, Lon} of a cell at level 24.
-spec cell_geometry(code()) -> [{float(), float()}].
cell_geometry(Code) -> cell_geometry(Code, ?MAX_LEVEL).

%% @doc Return the 6 corner coordinates {Lat, Lon} of a cell at the given level.
-spec cell_geometry(code(), level()) -> [{float(), float()}].
cell_geometry({Q, R}, Level) ->
    {CLat, CLon} = decode({Q, R}, Level),
    {CX, CY} = latlon_to_xy(CLat, CLon),
    BaseRadius = ?BQ_X / math:sqrt(3.0),
    Radius = BaseRadius * math:pow(2, ?MAX_LEVEL - Level),
    Angles = [30, 90, 150, 210, 270, 330],
    [xy_to_latlon(
        CX + Radius * math:cos(A * math:pi() / 180),
        CY + Radius * math:sin(A * math:pi() / 180)
     ) || A <- Angles].

display(Code) -> display(Code, ?MAX_LEVEL).

display({Q, R}, Level) ->
    Bits = interleave_bits(<<Q:24>>, <<R:24>>, <<>>),
    BitLen = Level * 2,
    <<Prefix:BitLen/bitstring, _/bitstring>> = Bits,
    PadLen = (8 - (BitLen rem 8)) rem 8,
    Padded = <<Prefix/bitstring, 0:PadLen>>,
    Hex = binary:encode_hex(Padded),
    CharLen = (BitLen + 3) div 4,
    binary:part(Hex, 0, CharLen).

parse(S) ->
    CharLen = byte_size(S),
    PaddedS = if CharLen rem 2 =:= 0 -> S; true -> <<S/binary, "0">> end,
    Bin = binary:decode_hex(PaddedS),
    LastChar = binary:last(S),
    IsOdd = (LastChar >= $0 andalso LastChar =< $3),
    BitLen = if IsOdd -> (CharLen * 4) - 2; true -> CharLen * 4 end,
    <<Bits:BitLen/bitstring, _/bitstring>> = Bin,
    {Qp, Rp} = deinterleave_bits(Bits, 0, 0),
    ActualLevel = BitLen div 2,
    Shift = ?MAX_LEVEL - ActualLevel,
    {Qp bsl Shift, Rp bsl Shift}.

%% ---------------------------------------------------------------------------
%% Bit Interleaving (internal)
%% ---------------------------------------------------------------------------

interleave_bits(<<Q:1, QRest/bitstring>>, <<R:1, RRest/bitstring>>, Acc) ->
    interleave_bits(QRest, RRest, <<Acc/bitstring, Q:1, R:1>>);
interleave_bits(<<>>, <<>>, Acc) -> Acc.

deinterleave_bits(<<Q:1, R:1, Rest/bitstring>>, QAcc, RAcc) ->
    deinterleave_bits(Rest, (QAcc bsl 1) bor Q, (RAcc bsl 1) bor R);
deinterleave_bits(<<>>, QAcc, RAcc) ->
    {QAcc, RAcc}.

%% ---------------------------------------------------------------------------
%% Geometry (internal)
%% ---------------------------------------------------------------------------

latlon_to_xy(Lat, Lon) ->
    CosLat = math:cos(Lat * math:pi() / 180.0),
    {(Lon - ?REF_LON) * ?M_PER_DEG_LAT * CosLat,
     (Lat - ?REF_LAT) * ?M_PER_DEG_LAT}.

xy_to_latlon(X, Y) ->
    Lat = ?REF_LAT + Y / ?M_PER_DEG_LAT,
    CosLat = math:cos(Lat * math:pi() / 180.0),
    {Lat, ?REF_LON + X / (?M_PER_DEG_LAT * CosLat)}.

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
