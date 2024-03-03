SET COUNT=0
:MyLoop
    IF "%COUNT%" == "40" GOTO EndLoop
    python .\genRandomFloat01.py
    SET /A COUNT+=1
    GOTO MyLoop
:EndLoop
