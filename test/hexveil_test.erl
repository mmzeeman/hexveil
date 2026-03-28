-module(hexveil_test).

-include_lib("eunit/include/eunit.hrl").

encode_decode_round_trip_test() ->
    Points = [
        {52.3026, 4.6889},  %% Amsterdam
        {-33.8688, 151.2093}, %% Sydney
        {40.7128, -74.0060}, %% New York
        {-23.5505, -46.6333}, %% São Paulo
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
        hexveil:encode(52.3026, 4.6889),
        hexveil:encode(0.0, 0.0)
    ],
    lists:foreach(
      fun(Code) ->
          Full = hexveil:display(Code),
          ?assertEqual(12, byte_size(Full)),
          ?assertEqual(Code, hexveil:parse(Full)),
          ?assertEqual(9, byte_size(hexveil:display(Code, 18))),
          ?assertEqual(9, byte_size(hexveil:display(Code, 17))),
          ?assertEqual(8, byte_size(hexveil:display(Code, 16))),
          ?assertEqual(1, byte_size(hexveil:display(Code, 1)))
      end,
      Codes).

display_prefix_property_test() ->
    Code = hexveil:encode(52.3026, 4.6889),
    %% In Base-16 (Hex), prefixes only hold for EVEN levels.
    S16 = hexveil:display(Code, 16),
    S18 = hexveil:display(Code, 18),
    S24 = hexveil:display(Code, 24),

    ?assert(binary_prefix(S16, S18)),
    ?assert(binary_prefix(S18, S24)).

coarsen_consistency_test() ->
    Code = hexveil:encode(52.3126, 4.6589),
    ?assertEqual(Code, hexveil:coarsen(Code, 24)),
    %% Use coarsen to find parent ID, and verify parse(display) returns the same ID.
    ?assertEqual(hexveil:coarsen(Code, 18), hexveil:coarsen(hexveil:parse(hexveil:display(Code, 18)), 18)),
    ?assertEqual(hexveil:coarsen(Code, 17), hexveil:coarsen(hexveil:parse(hexveil:display(Code, 17)), 17)),
    ?assertEqual(hexveil:coarsen(Code, 16), hexveil:coarsen(hexveil:parse(hexveil:display(Code, 16)), 16)).

are_nearby_test() ->
    Code = hexveil:encode(52.3026, 4.6889),
    SameCell = hexveil:encode(52.302601, 4.688901),
    DifferentButNearby = hexveil:encode(52.3035, 4.6895),

    ?assert(hexveil:are_nearby(Code, SameCell, 24)),
    ?assert(hexveil:are_nearby(Code, SameCell, 18)),
    ?assertNot(hexveil:are_nearby(Code, DifferentButNearby, 24)),
    ?assert(hexveil:are_nearby(Code, DifferentButNearby, 16)).

neighbors_test() ->
    Code = hexveil:encode(52.3026, 4.6889),
    Neighbors = hexveil:neighbors(Code),
    ?assertEqual(6, length(Neighbors)),
    ?assertEqual(
        lists:sort([
            {element(1, Code) + 1, element(2, Code)},
            {element(1, Code), element(2, Code) + 1},
            {element(1, Code) - 1, element(2, Code) + 1},
            {element(1, Code) - 1, element(2, Code)},
            {element(1, Code), element(2, Code) - 1},
            {element(1, Code) + 1, element(2, Code) - 1}
        ]),
        lists:sort(Neighbors)
    ).

cell_bounds_sanity_test() ->
    Code = hexveil:encode(52.3026, 4.6889),
    {MinLat, MinLon, MaxLat, MaxLon} = hexveil:cell_bounds(Code),
    {CLat, CLon} = hexveil:decode(Code),
    ?assert(MinLat < CLat),
    ?assert(MinLon < CLon),
    ?assert(MaxLat > CLat),
    ?assert(MaxLon > CLon),
    ?assert(MaxLat > MinLat),
    ?assert(MaxLon > MinLon).

binary_prefix(Prefix, Binary) ->
    Len = byte_size(Prefix),
    byte_size(Binary) >= Len andalso binary:part(Binary, 0, Len) =:= Prefix.
