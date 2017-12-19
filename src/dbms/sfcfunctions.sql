-- Reference
-- https://www.postgresql.org/docs/current/static/sql-createfunction.html

-- plpython examples
-- http://www.postgresonline.com/journal/index.php?/categories/51-plpython
-- http://www.postgresonline.com/journal/archives/100-PLPython-Part-2-Control-Flow-and-Returning-Sets.html


-- FIXME rename sfc to pysfc ???
-- FIXME annotate functions (stable, immutable, cost, ...)

-- use this to initialize the virtual environment inside the current session
-- SELECT workon('/pak2/usrdata/martijn/ais-env/'); -- << Place the code here inside site-packages
-- SELECT python_sys_path();

-- Hilbert encode / decode
--------------------------------------------------------------------------------
-- select hencode(array[1, 2]);
CREATE OR REPLACE FUNCTION hencode (arr integer[])
  RETURNS bigint
AS $$
    import sfc.hilbert
    return sfc.hilbert.encode(arr)
$$ LANGUAGE plpython2u;

-- select hdecode(10, 2);
CREATE OR REPLACE FUNCTION hdecode (val bigint, dims integer)
  RETURNS integer[]
AS $$
    import sfc.hilbert
    return sfc.hilbert.decode(val, dims)
$$ LANGUAGE plpython2u;


-- Morton encode / decode
--------------------------------------------------------------------------------
-- select nencode(array[1, 2]);
CREATE OR REPLACE FUNCTION nencode (arr integer[])
  RETURNS bigint
AS $$
    import sfc.morton_norder
    return sfc.morton_norder.encode(arr)
$$ LANGUAGE plpython2u;

-- select ndecode(6, 2);
CREATE OR REPLACE FUNCTION ndecode (val bigint, dims integer)
  RETURNS integer[]
AS $$
    import sfc.morton_norder
    return sfc.morton_norder.decode(val, dims)
$$ LANGUAGE plpython2u;


-- Range search
--------------------------------------------------------------------------------
-- function for Hilbert range search
DROP FUNCTION sfc_hquery;
CREATE OR REPLACE FUNCTION 
    sfc_hquery(lo integer[], hi integer[])
        RETURNS 
    TABLE(lower bigint, upper bigint)
AS $$
    import sfc.relate
    import sfc.query_hilbert
    return sfc.query_hilbert.hquery(query=sfc.relate.ndbox(lo, hi))
$$ LANGUAGE plpythonu;


-- function for N-order range search
DROP FUNCTION sfc_nquery;
CREATE OR REPLACE FUNCTION 
    sfc_nquery(lo integer[], hi integer[])
        RETURNS 
    TABLE(lower bigint, upper bigint)
AS $$
    import sfc.relate
    import sfc.query_norder
    return sfc.query_norder.nquery(query=sfc.relate.ndbox(lo, hi))
$$ LANGUAGE plpythonu;



-- Functions for scaling the values
--------------------------------------------------------------------------------
-- Scaling
DROP FUNCTION sfc_transform_scale;
DROP FUNCTION sfc_transform_dim;
DROP FUNCTION sfc_transform_dim_inv;
CREATE OR REPLACE FUNCTION sfc_transform_scale(old_lower numeric, old_upper numeric, new_lower numeric, new_upper numeric) RETURNS numeric AS $$
        DECLARE
            delta_old numeric;
            delta_new numeric;
        BEGIN
            delta_old := old_upper - old_lower;
            delta_new := new_upper - new_lower;
            RETURN delta_new / delta_old;
        END;
$$ LANGUAGE plpgsql
IMMUTABLE
;

-- select sfc_transform_scale(-180::numeric, 180::numeric, 0::numeric, pow(2,21)::numeric);
-- sfc_transform_scale 
-----------------------
--           5825.4224

CREATE OR REPLACE FUNCTION sfc_transform_dim(value numeric, translate numeric, scale numeric) RETURNS numeric AS $$
        BEGIN
            return (value - translate) * scale;
        END;
$$ LANGUAGE plpgsql
IMMUTABLE
;

CREATE OR REPLACE FUNCTION sfc_transform_dim_inv(value numeric, translate numeric, scale numeric) RETURNS numeric AS $$
        BEGIN
            return (value / scale) + translate;
        END;
$$ LANGUAGE plpgsql
IMMUTABLE
;

