@echo off
echo Establishing SSH tunnel...
start chrome http://localhost:8888

:: Prompt the user for the SSH username
set /p userInput="Enter your username: "
echo note: you will see no text when typing your password
:: Use the user input to replace "gebruiker" in the SSH command
ssh -L 8888:localhost:8787 -4 %userInput%@kgye0003l001.research.erasmusmc.nl