import scipy
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import pandas as pd
import glob
from scipy import interpolate

basis_data= pd.read_csv('G:/大量データ解析/データ格納庫/BLUE-JAY/BLUE-JAY_再探索結果_ver3.csv', header=None,engine='python')
basis=basis_data.values
"""サンプリング周波数"""
fs_HMS=33
fs_SIMS=10
fs_VDR=1
cut_filter=0.02

HMS_path='G:/大量データ解析/データ格納庫/BLUE-JAY/HMS_33Hz/*.csv'
SIMS_path='G:/大量データ解析/データ格納庫/BLUE-JAY/SIMS_10Hz/*.csv'
VDR_path='G:/大量データ解析/データ格納庫/BLUE-JAY/VDR/*.csv'

f1=open('G:/大量データ解析/データ格納庫/BLUE-JAY/BLUE-JAY_by_FFT.csv',"w")
f1.write("year,month,day,time,kn,HDG,H1,theta1,lamda1,w1,H2,theta2,lambda2,w2,error,HbyJWA,TbyJWA,keido,ido,P,R,H,S"+str('\n'))
f2=open('G:/大量データ解析/データ格納庫/BLUE-JAY/m0_BLUE-JAY_by_FFT.csv',"w")
f2.write("year,month,day,time,kn,HDG,H1,theta1,lamda1,w1,H2,theta2,lambda2,w2,error,HbyJWA,TbyJWA,keido,ido,P,R,H,S"+str('\n'))

HMS_list=[str(i[47:57]) for i in glob.glob(HMS_path)]
SIMS_list=[str(i[46:56]) for i in glob.glob(SIMS_path)]
VDR_list=[str(i[39:49]) for i in glob.glob(VDR_path)]
        
def Searching(date_number,listname):
    judge=False
    for i ,name in enumerate(listname):
        if name == date_number:
            judge=True
            break
    return i,judge
    
def HMS(Date):
    row, judge=Searching(Date,HMS_list)
    if judge==True:
        HMSdata=pd.read_csv(glob.glob(HMS_path)[row],engine='python',usecols=[5,6])
        HMSdata=HMSdata.values
        judge=False
        if len(HMSdata) >= fs_HMS*60*50:#当該csvでデータ量を確認
            judge=True
    return judge , np.array(HMSdata[:,0])/2+np.array(HMSdata[:,1])/2
    
def SIMS(Date):
    row,judge=Searching(Date,SIMS_list)
    if judge==True:
        SIMSdata=pd.read_csv(glob.glob(SIMS_path)[row],engine='python',usecols=[1,2,8])
        SIMSdata=SIMSdata.values
        judge=Flase
        if len(SIMSdata) >= fs_SIMS*60*50:#当該csvでデータ量を確認
            judge=True
    return judge , np.array(SIMSdata[:,0]), np.array(SIMSdata[:,1]), np.array(SIMSdata[:,2])
    
def VDR(Date):
    row , judge=Searching(Date,VDR_list)
    if judge==True:
        VDRdata=pd.read_csv(glob.glob(SIMS_path)[row],engine='python',usecols=[7,8])#status_water_speed,HDG
        VDRdata=VDRdata.values
        judge=False
        if len(VDRdata) >= fs_VDR*60*50:#当該csvでデータ量を確認
            if (np.max(VDRdata[:,0])-np.min(VDRdata[:,0])) <=2:#1時間で2kn以上の変化がないことを確認
                #角度変化の時系列作成
                direction_t=np.zeros(len(VDRdata))
                for i in range(len(VDRdata)-1):
                    direction_t[i+1] = np.arctan2(np.sin(VDRdata[i+1,1]-VDRdata[i,1]),np.cos(VDRdata[i+1,1]-VDRdata[i,1]))*180/np.pi+ direction_t[i]
                #1時間で船首方位が15°以内であることを確認
                if (np.max(direction_t)-np.min(direction_t)) <= 15 :
                    judge=True
        return judge , np.mean(VDRdata[:,0]), VDRdata[0,1]+np.mean(VDRdata[:,1])

def FFT(time_series,fs):
    freq=np.linspace(0,fs,len(time_series)
    dt=1/fs
    fn=1/dt/2
    F=np.fft.fft(time_series)
    F=F/(len(time_series)/2)
    F[(freq>=fn)]=0#ナイキスト周波数以降はカット
    F[(freq<=cut_filter)]=0#ロ―カット
    return F

def Integration(x,fs):#台形積分
    f=np.zeros(len(x))
    for i in range(len(x)-1):
        f[i+1]=(x[i+1]+x[i])*(1/fs)/2+f[i]
    return f

def Writing_basis(n_basis,Date,speed_mean,HDG_mean,f):
    #"year,month,day,time,kn,HDG,H1,theta1,lamda1,w1,H2,theta2,lambda2,w2,error,HbyJWA,TbyJWA,keido,ido
    f.write(str(Date[0:4])+","+str(Date[4:6])+","+str(Date[6:8])+","+str(Date[8:10]+","+str(speed_mean)+","+str(HDG_mean)+","+
        str(basis[n_basis][0])+","+str(basis[n_basis][1])+","+str(basis[n_basis][2])+","+str(basis[n_basis][3])+","+
        str(basis[n_basis][4])+","+str(basis[n_basis][5])+","+str(basis[n_basis][6])+","+str(basis[n_basis][7])+","+str(basis[n_basis][8])+","+
        str(basis[n_basis][17])+","+str(basis[n_basis][18])+","+str(basis[n_basis][23])+","+str(basis[n_basis][24])+","+)   

def Writing(x,f):
    for i in range(len(x)):
        f.write(str(x)+",")

def Writing_m0(x,f):
    m0=0
    for i in range(len(x)-1):
        m0+=(x[i]+x[i+1])*0.03*2
    f.write(str(m0)+",")
    
class Calculation:
    def __init__(self,pitch_t,roll_t,heave_t,GMPGMS_t):
        self.pitch_t =np.real(np.fft.ifft(FFT(pitch_t,fs_SIMS)))*len(pitch_t)
        self.roll_t  =np.real(np.fft.ifft(FFT(roll_t,fs_SIMS)))*len(roll_t)
        self.heave_t =self.Heave_calculation(heave_t)
        self.GMPGMS_t=np.real(np.fft.ifft(FFT(GMPGMS_t,fs_HMS)))*len(GMPGMS_t)
        self.correction_factor=1/0.5019
        self.encounter_omega=np.arange(0,2.11,0.03)
        
    def Heave_calculation(self,x):#Heave成分の時系列データ作成
        Z_acc=np.real(np.fft.ifft(FFT(x,fs_SIMS)))*len(x)#周波数カット
        Z_vel=Integration(Z_acc,fs_SIMS)#台形積分
        Z_vel=np.real(np.fft.ifft(FFT(Z_vel,fs_SIMS)))*len(Z_vel)#周波数カット
        Z_dis=Integration(Z_vel,fs_SIMS)#台形積分
        Z_dis=np.real(np.fft.ifft(FFT(Z_dis,fs_SIMS)))*len(Z_dis)#周波数カット
        heave_actual=Z_dis-49.55878*np.sin(self.pitch_t*np.pi/180)#ピッチ成分カット
        heave_actual=np.real(np.fft.ifft(FFT(heave_actual,fs_SIMS)))*len(heave_actual)#周波数カット
        return heave_actual

    def Framesize(self,number):
        return 2/301*number
    
    def Slide_size(self,number):
        return number/2
    
    def Hanning(self,number):
        return np.array([0.5-0.5*np.cos(2*np.pi/(self.Framesize(number)-1)*i) for i in range(self.Framesize(number))])
        
    def Interpolate(self,number,PSD1):
        S=interpolate.interp1d(freq_frame[0:int(self.Framesize(number)/2)+1]*2*np.pi, PSD1,kind='nearset')
        return S(self.encounter_omega)
    
    def Welch_ModifiedPeriodogram(self,X,Y):
        PSD_all=np.zeros(len(self.encounter_omega))
        for i in range(300):
            X_taper=X[i*self.Slide_size(len(X)) : i*self.Slide_size(len(X))+self.Framesize(len(X))]*self.Hanning(len(X))
            Y_taper=Y[i*self.Slide_size(len(Y)) : i*self.Slide_size(len(Y))+self.Framesize(len(Y))]*self.Hanning(len(Y))
            """FFT実行"""
            F_X=np.fft.fft(X_taper)*self.correction_factor
            F_Y=np.fft.fft(Y_taper)*self.correction_factor
            """内挿"""
            F_X=self.Interpolate(len(X),F_X)
            F_Y=self.Interpolate(len(Y),F_Y)
            """掛け算"""
            PSD=F_X.conjugate()*F_Y
            PSD[1:int(self.Framesize(len(self.encounter_omega))/2)-1]=2*PSD[1:int(self.Framesize(len(self.encounter_omega))/2)-1].copy()
            PSD_all += PSD
        return PSD_all/300
  
def main():
    for n_basis in range(len(basis)):
        Date=str(basis[n_basis][9])+str(basis[n_basis][10])+str('%02.0f'%(basis[n_basis][11]))+str('%02.0f'%(basis[n_basis][12]))
        
        VDRjudge, spped_mean , HDG_mean = VDR(Date)#VDRデータから対象csvファイル検索
        if VDRjudge == False:
            continue
        
        SIMSjudge, pitch_t, roll_t ,heave_t = SIMS(Date)#SIMSデータから対象csvファイル検索
        if SIMSjudge == False:
            continue
        
        HMSjudge, GMPGMS_t = csv_judgement.HMS(Date)#HMSデータから対象csvファイル検索
        if HMSjudge==False:
            continue
        
        DataSet=Calculation(pitch_t,roll_t,heave_t,GMPGMS_t)
        """計算&書き込み開始"""
        Writing_basis(n_basis,Date,speed_mean,HDG_mean,f1)
        Writing_basis(n_basis,Date,speed_mean,HDG_mean,f2)
        XY = DataSet.Welch_ModifiedPeriodogram(pitch_t,pitch_t)
        Writing(np.abs(XY),f1)
        Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(roll_t,roll_t)
        Writing(np.abs(XY),f1)
        Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(heave_t,heave_t)
        Writing(np.abs(XY),f1)
        Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,GMPGMS_t)
        Writing(np.abs(XY),f1)
        Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(pitch_t,roll_t)
        Writing(XY.real,f1)
        Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(roll_t,heave_t)
        Writing(XY.real,f1)
        Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(heave_t,pitch_t)
        Writing(XY.real,f1)
        Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,roll_t)
        Writing(XY.real,f1)
        Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,heave_t)
        Writing(XY.real,f1)
        Writing(XY.imag,f1)        
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,pitch_t)
        Writing(XY.real,f1)
        Writing(XY.imag,f1)
        
        f1.write('/n')
        f2.write('/n')
        print(n_basis)
    print("End")    

if __name__ == '__main__':
    main()
