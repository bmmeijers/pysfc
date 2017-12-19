
-- select floor(min()) from kdz_ais_bits;
-- select ceil(max()) from kdz_ais_bits;

-- SELECT workon('/pak2/usrdata/martijn/ais-env/'); -- << Place the code here inside site-packages
-- SELECT python_sys_path();


create table kdz_1day as
with transform_params as(
    select
    sfc_transform_scale(-180, 180, 0, pow(2,21)::numeric) as scale_east, 
    sfc_transform_scale( -90,  90, 0, pow(2,21)::numeric) as scale_north,
    sfc_transform_scale(1475272799, 1483225200, 0, pow(2,21)::numeric) as scale_time,
    -180 as translate_east,
     -90 as translate_north,
    1475272799.0 as translate_time
)
select hencode(array[scaled_east, scaled_north, scaled_time]) as hilbert, ts, pt, pyl as payload
from
(
select
    sfc_transform_dim(ais_easting(pyl)::numeric, transform_params.translate_east, transform_params.scale_east)::int as scaled_east,
    sfc_transform_dim(ais_northing(pyl)::numeric, transform_params.translate_north, transform_params.scale_north)::int as scaled_north,
    sfc_transform_dim(extract(epoch from ts)::numeric, transform_params.translate_time, transform_params.scale_time)::int as scaled_time,
    ts,
    ais_point(pyl) as pt,
    pyl
from (select ts, payload as pyl from kdz_ais_bits where ts < '2016-10-01 23:59:59' order by ts) as subset left join transform_params on TRUE
) as scaled_data
;



-- Prepare table by subsetting the data
create table kdz_1day as
with ais_data as(
    select
    sfc_transform_scale(-180, 180, 0, pow(2,21)::numeric) as scale_east, 
    sfc_transform_scale( -90,  90, 0, pow(2,21)::numeric) as scale_north,
    sfc_transform_scale(1475272799, 1483225200, 0, pow(2,21)::numeric) as scale_time,
    -180 as translate_east,
     -90 as translate_north,
    1475272799.0 as translate_time
)
select
    hencode(array[scaled_east, scaled_north, scaled_time]) as hilbert, ts, pt, pyl as payload
from
(
    select
        sfc_transform_dim(ais_easting(pyl)::numeric, ais_data.translate_east, ais_data.scale_east)::int as scaled_east,
        sfc_transform_dim(ais_northing(pyl)::numeric, ais_data.translate_north, ais_data.scale_north)::int as scaled_north,
        sfc_transform_dim(extract(epoch from ts)::numeric, ais_data.translate_time, ais_data.scale_time)::int as scaled_time,
        ts,
        ais_point(pyl) as pt,
        pyl
    from 
    (
        select ts, payload as pyl from kdz_ais_bits where ts < '2016-10-01 23:59:59Z' order by ts
    ) as subset
    left join
        ais_data on TRUE
) as scaled_data
;

-- make the index
CREATE INDEX
    i__kdz_1day__hilbert
ON
    kdz_1day (hilbert asc)
TABLESPACE
    indx
;

-- cluster on the index
CLUSTER kdz_1day USING i__kdz_1day__hilbert
;

-- vacuum + analyze the table
VACUUM ANALYZE kdz_1day;



-- POLYGON((3 51, 3 54, 6 54, 6 51, 3 51))

-- make range table

drop table kdz_1day_ranges;

-- create unlogged
create temp table kdz_1day_ranges as select * from (
with transform_params as(
    select
        sfc_transform_scale(-180, 180, 0, pow(2,21)::numeric) as scale_east, 
        sfc_transform_scale( -90,  90, 0, pow(2,21)::numeric) as scale_north,
        sfc_transform_scale(1475272799, 1483225200, 0, pow(2,21)::numeric) as scale_time,
        -180 as translate_east,
        -90 as translate_north,
        1475272799.0 as translate_time
    )
select * from sfc_hquery(
    array[
            sfc_transform_dim(3,          (select translate_east  from transform_params), (select scale_east  from transform_params))::int,
            sfc_transform_dim(51,         (select translate_north from transform_params), (select scale_north from transform_params))::int,
            sfc_transform_dim(1475280000, (select translate_time  from transform_params), (select scale_time  from transform_params))::int
    ],
    array[
            sfc_transform_dim(6,          (select translate_east  from transform_params), (select scale_east  from transform_params))::int,
            sfc_transform_dim(54,         (select translate_north from transform_params), (select scale_north from transform_params))::int,
            sfc_transform_dim(1475280600, (select translate_time  from transform_params), (select scale_time  from transform_params))::int
    ]
) as foo)
as foofoo;


create index i__kdz_1day_ranges__lower on kdz_1day_ranges (lower asc);
cluster kdz_1day_ranges using i__kdz_1day_ranges__lower;
vacuum analyze kdz_1day_ranges;

explain analyze select * from kdz_1day d, kdz_1day_ranges r where d.hilbert >= r.lower and d.hilbert < r.upper;

drop table if exists kdz_1day_subset;
create table kdz_1day_subset as select * from kdz_1day d, kdz_1day_ranges r where d.hilbert >= r.lower and d.hilbert < r.upper;

