-module(isea4h_tests).
-include_lib("eunit/include/eunit.hrl").

%% Round-trip encode/decode at origin
roundtrip_test() ->
    Code = isea4h:encode(0.0, 0.0, 7),
    {Lon, Lat} = isea4h:decode(Code),
    ?assert(math:abs(Lon - 0.0) < 0.01),
    ?assert(math:abs(Lat - 0.0) < 0.01).

%% encode/2 should default to resolution 7
default_res_test() ->
    Code = isea4h:encode(1.23, 4.56),
    [_,Digits] = string:split(binary_to_list(Code), "-"),
    ?assertEqual(7, length(Digits)).

%% parent should remove one digit at the end (or leave as-is for resolution 1)
parent_test() ->
    Code = isea4h:encode(10.0, 20.0, 6),
    Parent = isea4h:parent(Code),
    [_,Digits] = string:split(binary_to_list(Code), "-"),
    [_,Pdigits] = string:split(binary_to_list(Parent), "-"),
    ?assertEqual(length(Digits)-1, length(Pdigits)).

%% neighbors returns six binary neighbor codes
neighbors_test() ->
    Code = isea4h:encode(10.0, 20.0, 5),
    N = isea4h:neighbors(Code),
    ?assertEqual(6, length(N)),
    lists:foreach(fun(C) -> ?assert(is_binary(C)) end, N).

%% privacy_code returns codes of the expected length for each privacy level
privacy_test() ->
    [_,D1] = string:split(binary_to_list(isea4h:privacy_code(0.0,0.0, city)), "-"),
    [_,D2] = string:split(binary_to_list(isea4h:privacy_code(0.0,0.0, district)), "-"),
    [_,D3] = string:split(binary_to_list(isea4h:privacy_code(0.0,0.0, neighbourhood)), "-"),
    [_,D4] = string:split(binary_to_list(isea4h:privacy_code(0.0,0.0, block)), "-"),
    ?assertEqual(4, length(D1)),
    ?assertEqual(6, length(D2)),
    ?assertEqual(7, length(D3)),
    ?assertEqual(9, length(D4)).

%% ico_verts returns the expected number of vertices (12)
ico_test() ->
    V = isea4h:ico_verts(),
    ?assertEqual(12, length(V)).
