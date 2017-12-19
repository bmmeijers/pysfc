-- Functions to activate a virtual environment for the Python interpreter

-- https://gist.github.com/Gwill/071e4c475dd013a97a2e
CREATE OR REPLACE FUNCTION workon(venv text)
  RETURNS void AS
$BODY$
    import os
    import sys
    
    if sys.platform in ('win32', 'win64', 'cygwin'):
        activate_this = os.path.join(venv, 'Scripts', 'activate_this.py')
    else:
        activate_this = os.path.join(venv, 'bin', 'activate_this.py')
        
    exec(open(activate_this).read(), dict(__file__=activate_this))
$BODY$
LANGUAGE plpython2u VOLATILE;


--
CREATE OR REPLACE FUNCTION python_sys_path()
  RETURNS text[] AS
$BODY$
    import sys
    return sys.path
$BODY$
  LANGUAGE plpython2u VOLATILE;
