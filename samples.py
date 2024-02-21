DOG_SAMPLES = '''
SAMN06172495
SAMN06172481
SAMN06172523
SAMN06172473
SAMN06172512
SAMN06172522
SAMN06172491
SAMN06172421
SAMN06172511
SAMN06172505
SAMN06172453
SAMN06172443
'''

DOG_SAMPLES = [ell.strip() for ell in DOG_SAMPLES.strip().split('\n')]
DOG_STUDY = 'Coelho_dogs_2018'

ZELLER_SAMPLES = '''
SAMEA4817960
SAMEA4817939
SAMEA4817907
SAMEA4817900
SAMEA4817972
SAMEA4817969
SAMEA4817928
SAMEA4817953
SAMEA4817943
SAMEA4817947
'''

ZELLER_SAMPLES = [('PAIRED/'+ell.strip()) for ell in ZELLER_SAMPLES.strip().split('\n')]
ZELLER_STUDY = 'Zeller_CRC_2014'

