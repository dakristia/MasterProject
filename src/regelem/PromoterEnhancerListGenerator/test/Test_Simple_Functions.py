import sys
import pytest #type:ignore
sys.path.insert(0, '..')
from Functions import *


def test_check_in_range():

    #Both numbers are in range of atleast one range set
    #Should all be True
    assert check_in_range(3000,3300,0,5000,5000,10000) == True
    assert check_in_range(3000,3300,0,5000,-1,-1) == True
    assert check_in_range(3000,3300,-1,-1,0,5000) == True

    #Should all be False
    assert check_in_range(3000,3300,-1,-1,-1,-1) == False
    assert check_in_range(3000,3300,5000,10000,10000,15000) == False

    #Only one number is in range of atleast one range set
    assert check_in_range(3000,3300,2999,3200,-1,-1) == True
    assert check_in_range(3000,3300,-1,-1,2999,3200) == True
    assert check_in_range(3000,3300,-1,-1,3200,3400) == True

    #Checking edge cases
    #Should all be true
    assert check_in_range(3000,3300,3000,3010,-1,-1)
    assert check_in_range(3000,3300,3300,3310,-1,-1)
    assert check_in_range(3000,3300,2500,3000,-1,-1)
    assert check_in_range(3000,3300,3200,3300,-1,-1)

    assert check_in_range(3000,3300,3301,3330,-1,-1) == False
    assert check_in_range(3000,3300,2000,2999,-1,-1) == False

def test_nums_to_str_range():

    assert nums_to_str_range(1000,6000) == "1000-6000"

def test_get_digit():
    assert get_digit(54321,0) == 1
    assert get_digit(54321,1) == 2
    assert get_digit(54321,2) == 3
    assert get_digit(54321,3) == 4
    assert get_digit(54321,4) == 5
    assert get_digit(54321,5) == 0


def test_round_up_and_down():
    down,up = round_up_and_down(3000,5000)
    assert down == 0
    assert up == 5000

    down,up = round_up_and_down(3000,10000)
    assert down == 0
    assert up == 10000

    down,up = round_up_and_down(3000,30)
    assert down == 3000
    assert up == 3030

    down,up = round_up_and_down(3100,200)
    assert down == 3000
    assert up == 3200
