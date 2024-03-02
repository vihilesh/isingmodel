SET COUNT=0
:MyLoop
    IF "%COUNT%" == "40" GOTO EndLoop
    python .\ANUQRNG.py
    SET /A COUNT+=1
    GOTO MyLoop
:EndLoop
