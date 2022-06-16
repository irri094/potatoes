@echo off
if exist log.out del log.out
echo starting
for /l %%x in (1, 1, %1) do (
    TestGenerator > input.in
    solution < input.in > output.out 2>> log.out
    brute < input.in > output2.out 2>> log.out
    fc output.out output2.out > diagnostics || exit /b
    echo %%x
)
echo all tests passed
