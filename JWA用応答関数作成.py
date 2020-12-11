import pandas as pd
import numpy as np
from scipy import interpolate

strip=pd.read_csv("20knot_wide_ver2.csv",engine='python',header=None,skiprows=1,
usecols=[0,3,9,11,13,17])
strip=strip.values
print(strip[0])
print(strip.shape)

"""諸々の定義"""
n_lambda  =80
n_angle   =72
omegalist =np.array(strip[:n_lambda,1])
freq_JWA=2*np.pi*np.array([0.0445953,0.0486315,0.053033,0.0578329,0.0630672,0.0687753,0.075,0.0817881,
               0.0891906,0.097263,0.106066,0.1156658,0.1261345,0.1375506,0.15,0.1635762,0.1783811,
               0.194526,0.2121321,0.2313317,0.252269,0.2751013,0.3000001,0.3271524,0.3567623])
freq_JWA=freq_JWA[::-1]

"""csvファイルをオープン"""
f=open("strip_on_JWA.csv","w")
f.write("Fn,freq,heave_amp,roll_amp,pitch_amp,stress_amp")
f.write('\n')

"""for分を回してDataFrameを作成"""
all_data=[[],[],[],[]]
for i in range(n_angle):
    """内挿条件の作成
    F_pitch =interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,4])
    F_roll  =interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,3])
    F_heave =interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,2])
    F_stress=interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,5])"""
    
    """内挿を実行してデータフレームに格納"""
    all_data[0]=all_data[0]+interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,4])(freq_JWA)
    all_data[1]=all_data[1]+interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,3])(freq_JWA)
    all_data[2]=all_data[2]+interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,2])(freq_JWA)
    all_data[3]=all_data[3]+interpolate.inter1d(omegalist,strip[i*n_lambda:(i+1)*n_lambda,5])(freq_JWA)

print(all_data.T.shape)

"""書き込み"""
for i in range(n_lambda*n_angle):
    f.write(str(freq_JWA[int(i % n_lambda)])+",")
    for j in all_data.T[i]:
        f.write(str(j)+",")
    f.write('\n')
f.close()
print("end")
