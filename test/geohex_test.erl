-module(geohex_test).

-include_lib("eunit/include/eunit.hrl").

encode_decode_round_trip_test() ->
    Points = [
        {52.3026,   4.6889},   %% Amsterdam
        {52.0907,   5.1214},   %% Utrecht
        {51.9244,   4.4777},   %% Rotterdam
        {52.3676,   4.9041},
        {51.4416,   5.4697},
        {53.2194,   6.5665},
        {40.7128,  -74.0060},  %% New York
        {35.6762,  139.6503},  %% Tokyo
        {-33.8688, 151.2093},  %% Sydney
        {-23.5505, -46.6333}   %% São Paulo
    ],
    lists:foreach(
      fun({Lat, Lon}) ->
          Code = geohex:encode(Lat, Lon),
          {Lat2, Lon2} = geohex:decode(Code),
          ?assert(is_integer(element(1, Code))),
          ?assert(is_integer(element(2, Code))),
          ?assert(abs(Lat - Lat2) < 0.02),
          ?assert(abs(Lon - Lon2) < 0.02)
      end,
      Points),
    %% Edge cases: encode/decode must not crash (projection distortion expected)
    EdgeCases = [
        {-90.0,   0.0},   %% South Pole
        { 90.0,   0.0},   %% North Pole
        {  0.0,  180.0},  %% anti-meridian east
        {  0.0, -180.0}   %% anti-meridian west
    ],
    lists:foreach(
      fun({Lat, Lon}) ->
          Code = geohex:encode(Lat, Lon),
          ?assert(is_integer(element(1, Code))),
          ?assert(is_integer(element(2, Code))),
          _ = geohex:decode(Code)
      end,
      EdgeCases).

display_parse_round_trip_test() ->
    Codes = [
        geohex:encode(52.3026, 4.6889),
        geohex:encode(52.0907, 5.1214),
        geohex:encode(51.9244, 4.4777)
    ],
    lists:foreach(
      fun(Code) ->
          Full = geohex:display(Code),
          ?assertEqual(10, byte_size(Full)),
          ?assertEqual(Code, geohex:parse(Full)),
          ?assertEqual(8, byte_size(geohex:display(Code, 4))),
          ?assertEqual(6, byte_size(geohex:display(Code, 3))),
          ?assertEqual(4, byte_size(geohex:display(Code, 2))),
          ?assertEqual(2, byte_size(geohex:display(Code, 1)))
      end,
      Codes).

display_prefix_property_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    S1 = geohex:display(Code, 1),
    S2 = geohex:display(Code, 2),
    S3 = geohex:display(Code, 3),
    S4 = geohex:display(Code, 4),
    S5 = geohex:display(Code, 5),
    ?assertEqual(true, binary_prefix(S1, S2)),
    ?assertEqual(true, binary_prefix(S2, S3)),
    ?assertEqual(true, binary_prefix(S3, S4)),
    ?assertEqual(true, binary_prefix(S4, S5)).

coarsen_consistency_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    ?assertEqual(Code, geohex:coarsen(Code, 5)),
    ?assertEqual(geohex:coarsen(Code, 4), geohex:parse(geohex:display(Code, 4))),
    ?assertEqual(geohex:coarsen(Code, 3), geohex:parse(geohex:display(Code, 3))),
    ?assertEqual(geohex:coarsen(Code, 2), geohex:parse(geohex:display(Code, 2))),
    ?assertEqual(geohex:coarsen(Code, 1), geohex:parse(geohex:display(Code, 1))).

are_nearby_test() ->
    Code = geohex:encode(52.3026, 4.6889),
    SameCell = geohex:encode(52.3027, 4.6890),
    %% ~3km away — guaranteed different level-5 cell at 2.5km resolution
    FarAway = geohex:encode(52.3300, 4.7200),

    ?assert(geohex:are_nearby(Code, SameCell, 5)),
    ?assert(geohex:are_nearby(Code, SameCell, 4)),
    ?assert(geohex:are_nearby(Code, SameCell, 3)),
    ?assertNot(geohex:are_nearby(Code, FarAway, 5)),
    ?assert(geohex:are_nearby(Code, FarAway, 4) orelse
             geohex:are_nearby(Code, FarAway, 3)).

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
