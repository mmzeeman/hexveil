-module(hexveil_test).

-include_lib("eunit/include/eunit.hrl").

encode_decode_round_trip_test() ->
    Points = [
        {52.3616, 4.8784},  %% Vondelpark
        {-33.8688, 151.2093}, %% Sydney
        {40.7128, -74.0060}, %% New York
        {0.0, 0.0}            %% Null Island
    ],
    lists:foreach(
      fun({Lat, Lon}) ->
          Code = hexveil:encode(Lat, Lon),
          {Lat2, Lon2} = hexveil:decode(Code),
          ?assert(is_integer(element(1, Code))),
          ?assert(is_integer(element(2, Code))),
          ?assert(abs(Lat - Lat2) < 0.0001),
          ?assert(abs(Lon - Lon2) < 0.0001)
      end,
      Points).

display_parse_round_trip_test() ->
    Codes = [
        hexveil:encode(52.3616, 4.8784),
        hexveil:encode(0.0, 0.0)
    ],
    lists:foreach(
      fun(Code) ->
          Full = hexveil:display(Code),
          ?assertEqual(12, byte_size(Full)),
          ?assertEqual(Code, hexveil:parse(Full)),
          ?assertEqual(9, byte_size(hexveil:display(Code, 18))),
          ?assertEqual(8, byte_size(hexveil:display(Code, 16)))
      end,
      Codes).

display_prefix_property_test() ->
    Code = hexveil:encode(52.3616, 4.8784),
    S16 = hexveil:display(Code, 16),
    S18 = hexveil:display(Code, 18),
    S24 = hexveil:display(Code, 24),
    ?assert(binary_prefix(S16, S18)),
    ?assert(binary_prefix(S18, S24)).

coarsen_consistency_test() ->
    Code = hexveil:encode(52.3616, 4.8784),
    ?assertEqual(Code, hexveil:coarsen(Code, 24)),
    ?assertEqual(hexveil:coarsen(Code, 18), hexveil:coarsen(hexveil:parse(hexveil:display(Code, 18)), 18)),
    ?assertEqual(hexveil:coarsen(Code, 16), hexveil:coarsen(hexveil:parse(hexveil:display(Code, 16)), 16)).

are_nearby_test() ->
    Code = hexveil:encode(52.3616, 4.8784),
    SameCell = hexveil:encode(52.361601, 4.878401),
    DifferentButNearby = hexveil:encode(52.3650, 4.8850),
    ?assert(hexveil:are_nearby(Code, SameCell, 24)),
    ?assert(hexveil:are_nearby(Code, SameCell, 18)),
    ?assertNot(hexveil:are_nearby(Code, DifferentButNearby, 24)),
    ?assert(hexveil:are_nearby(Code, DifferentButNearby, 14)).

binary_prefix(Prefix, Binary) ->
    Len = byte_size(Prefix),
    byte_size(Binary) >= Len andalso binary:part(Binary, 0, Len) =:= Prefix.
