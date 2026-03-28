-module(geohex_test).

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
          Code = geohex:encode(Lat, Lon),
          {Lat2, Lon2} = geohex:decode(Code),
          ?assert(is_integer(element(1, Code))),
          ?assert(is_integer(element(2, Code))),
          ?assert(abs(Lat - Lat2) < 0.01),
          ?assert(abs(Lon - Lon2) < 0.01)
      end,
      Points).

display_parse_round_trip_test() ->
    Codes = [
        geohex:encode(52.3026, 4.6889),
        geohex:encode(0.0, 0.0)
    ],
    lists:foreach(
      fun(Code) ->
          Full = geohex:display(Code),
          ?assertEqual(7, byte_size(Full)),
          ?assertEqual(Code, geohex:parse(Full)),
          ?assertEqual(6, byte_size(geohex:display(Code, 6))),
          ?assertEqual(5, byte_size(geohex:display(Code, 5))),
          ?assertEqual(4, byte_size(geohex:display(Code, 4))),
          ?assertEqual(1, byte_size(geohex:display(Code, 1)))
      end,
      Codes).

display_prefix_property_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    S1 = geohex:display(Code, 1),
    S4 = geohex:display(Code, 4),
    S6 = geohex:display(Code, 6),
    S7 = geohex:display(Code, 7),

    %% In Base-49, each character represents one level exactly.
    ?assert(binary_prefix(S1, S4)),
    ?assert(binary_prefix(S4, S6)),
    ?assert(binary_prefix(S6, S7)).

coarsen_consistency_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    ?assertEqual(Code, geohex:coarsen(Code, 7)),
    %% parse(display(Code, L)) returns the level-7 coordinate of the cell's start.
    ?assertEqual(geohex:coarsen(Code, 6), geohex:coarsen(geohex:parse(geohex:display(Code, 6)), 6)),
    ?assertEqual(geohex:coarsen(Code, 4), geohex:coarsen(geohex:parse(geohex:display(Code, 4)), 4)),
    ?assertEqual(geohex:coarsen(Code, 1), geohex:coarsen(geohex:parse(geohex:display(Code, 1)), 1)).

are_nearby_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    SameCell = geohex:encode(52.3027, 4.6890),
    DifferentButNearby = geohex:encode(52.3040, 4.6910),

    ?assert(geohex:are_nearby(Code, SameCell, 7)),
    ?assert(geohex:are_nearby(Code, SameCell, 6)),
    ?assertNot(geohex:are_nearby(Code, DifferentButNearby, 7)),
    %% At level 4 (city scale), they should definitely be in the same cell.
    ?assert(geohex:are_nearby(Code, DifferentButNearby, 4)).

neighbors_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    Neighbors = geohex:neighbors(Code),
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
    Code = geohex:encode(52.3026, 4.6889),
    {MinLat, MinLon, MaxLat, MaxLon} = geohex:cell_bounds(Code),
    {CLat, CLon} = geohex:decode(Code),
    ?assert(MinLat < CLat),
    ?assert(MinLon < CLon),
    ?assert(MaxLat > CLat),
    ?assert(MaxLon > CLon),
    ?assert(MaxLat > MinLat),
    ?assert(MaxLon > MinLon).

binary_prefix(Prefix, Binary) ->
    Len = byte_size(Prefix),
    byte_size(Binary) >= Len andalso binary:part(Binary, 0, Len) =:= Prefix.
