-module(isea4h_tests).
-include_lib("eunit/include/eunit.hrl").

%% Round-trip encode/decode at origin
roundtrip_test() ->
    Locations = [
        {0.0, 0.0},
        {10.0, 20.0},
        {-30.0, 45.0},
        {120.0, -10.0},
        {0.0, 90.0}, % North Pole
        {0.0, -90.0} % South Pole
    ],
    lists:foreach(fun({Lon, Lat}) ->
        Res = 7,
        Code = isea4h:encode(Lon, Lat, Res),
        {DLon, DLat} = isea4h:decode(Code),
        MaxErr = 1.0,
        IsPole = abs(Lat) > 89.0,
        LonMatch = IsPole orelse abs(DLon - Lon) < MaxErr orelse abs(abs(DLon - Lon) - 360.0) < MaxErr,
        LatMatch = abs(DLat - Lat) < MaxErr,
        case LonMatch andalso LatMatch of
            true -> ok;
            false -> 
                io:format(user, "~nAt (~p, ~p) got (~p, ~p) code ~p~n", [Lon, Lat, DLon, DLat, Code]),
                ?assert(false)
        end
    end, Locations).

roundtrip_res_test() ->
    Lon = 10.0, Lat = 20.0,
    lists:foreach(fun(Res) ->
        Code = isea4h:encode(Lon, Lat, Res),
        {DLon, DLat} = isea4h:decode(Code),
        %% At Res 1, error can be large (half a face).
        MaxErr = 2.0 / math:pow(2.0, Res-5), %% Heuristic
        ActualMaxErr = lists:max([1.0, MaxErr]),
        case abs(DLon - Lon) < ActualMaxErr andalso abs(DLat - Lat) < ActualMaxErr of
            true -> ok;
            false -> 
                io:format(user, "~nAt res ~p got (~p, ~p) code ~p (MaxErr ~p)~n", [Res, DLon, DLat, Code, ActualMaxErr]),
                ?assert(false)
        end
    end, lists:seq(1, 12)).

%% encode/2 should default to resolution 7
default_res_test() ->
    Code = isea4h:encode(1.23, 4.56),
    [_,Digits] = string:split(binary_to_list(Code), "-"),
    ?assertEqual(7, length(Digits)).

%% parent should remove one digit at the end (or leave as-is for resolution 1)
parent_test() ->
    Code = isea4h:encode(10.0, 20.0, 6),
    Parent = isea4h:parent(Code),
    [_, Digits] = string:split(binary_to_list(Code), "-"),
    [_, Pdigits] = string:split(binary_to_list(Parent), "-"),
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

%% Verify the specific Lon/Lat of the icosahedron vertices
ico_coords_test() ->
    %% We need to reach into the internal to_xyz/decode logic or just verify the Verts list.
    %% Since ico_verts() returns XYZ, let's convert them back to Lon/Lat for verification.
    D2R = math:pi() / 180.0,
    Verts = isea4h:ico_verts(),
    Coords = [begin
                R = math:sqrt(X*X + Y*Y + Z*Z),
                Lat = math:asin(Z/R) / D2R,
                Lon = math:atan2(Y, X) / D2R,
                {Lon, Lat}
              end || {X, Y, Z} <- Verts],
    
    %% Expected Structure:
    %% 0: North Pole (90)
    %% 1-5: Upper Ring (~26.56 lat, 0, 72, 144, 216, 288 lon)
    %% 6-10: Lower Ring (~ -26.56 lat, 36, 108, 180, 252, 324 lon)
    %% 11: South Pole (-90)
    
    {_, NP_Lat} = lists:nth(1, Coords),
    ?assert(abs(NP_Lat - 90.0) < 0.0001),
    
    UpperRing = lists:sublist(Coords, 2, 5),
    lists:foreach(fun({_Lon, Lat}) ->
        ?assert(abs(Lat - 26.56505) < 0.001)
    end, UpperRing),
    
    LowerRing = lists:sublist(Coords, 7, 5),
    lists:foreach(fun({_Lon, Lat}) ->
        ?assert(abs(Lat + 26.56505) < 0.001)
    end, LowerRing),
    
    {_, SP_Lat} = lists:nth(12, Coords),
    ?assert(abs(SP_Lat + 90.0) < 0.0001),
    
    %% Check longitudes of upper ring (0, 72, 144, 216, 288)
    ExpectedUpperLon = [0.0, 72.0, 144.0, -144.0, -72.0], %% atan2 returns -180..180
    ActualUpperLon = [L || {L, _} <- UpperRing],
    lists:zipwith(fun(E, A) -> ?assert(abs(E - A) < 0.0001) end, ExpectedUpperLon, ActualUpperLon),
    
    ok.
