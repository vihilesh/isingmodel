SET COUNT=0
:MyLoop
    IF "%COUNT%" == "2000" GOTO EndLoop
    python .\ANUQRNG.py
    SET /A COUNT+=1
    timeout /t 2 /nobreak > NUL
    GOTO MyLoop
:EndLoop
./