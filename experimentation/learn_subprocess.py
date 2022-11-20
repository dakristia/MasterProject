import sys
import subprocess

def test_process():
    text1 = "hallo" 
    text2 = "der"

    with open("testfil.txt","w") as f:
       f.write("aa" + " " + "bb")

    return text1 + " " + text2


print(test_process())

print("allo der")
