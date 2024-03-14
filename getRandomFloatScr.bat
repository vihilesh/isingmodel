SET COUNT=0
:MyLoop
    IF "%COUNT%" == "5000" GOTO EndLoop
    python .\genRandomFloat01.py
    SET /A COUNT+=1
    timeout /t 1 /nobreak > NUL
    GOTO MyLoop
:EndLoop
