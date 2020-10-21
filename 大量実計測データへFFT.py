import scipy
import numpy as np
from scipy import signal
import pandas as pd
import glob
from scipy import interpolate

basis_data= pd.read_csv('G:/大量データ解析/データ格納庫/BLUE-JAY/BLUE-JAY_再探索結果_ver3.csv', header=None,engine='python')
basis=basis_data.values
print(np.shape(basis))
print(str(basis[0][9])+str(basis[0][10])+str('%02.0f'%(basis[0][11]))+str('%02.0f'%(basis[0][12])))

"""サンプリング周波数"""
fs_HMS=33
fs_SIMS=10
fs_VDR=1
cut_filter=0.02 #カット周波数

"""各データセットのpath"""
HMS_path='G:/大量データ解析/データ格納庫/BLUE-JAY/HMS_33Hz/*.csv'
SIMS_path='G:/大量データ解析/データ格納庫/BLUE-JAY/SIMS_10Hz/*.csv'
VDR_path='G:/大量データ解析/データ格納庫/BLUE-JAY/VDR/*.csv'
glob.glob(HMS_path)[0][47:57]

"""各データの日付情報リスト作成"""
HMS_list=[str(i[47:57]) for i in glob.glob(HMS_path)]
print(HMS_list[0])
SIMS_list=[str(i[46:56]) for i in glob.glob(SIMS_path)]
print(SIMS_list[0])
VDR_list=[str(i[39:49]) for i in glob.glob(VDR_path)]
print(VDR_list[0])

"""書き込みcsvファイル作成"""
f1=open('G:/大量データ解析/データ格納庫/BLUE-JAY/BLUE-JAY_by_FFT.csv',"w")
f1.write("year,month,day,time,kn,HDG,H1,theta1,lamda1,w1,H2,theta2,lambda2,w2,error,HbyJWA,TbyJWA,keido,ido,P,R,H,S,")
f1.write('\n')
f2=open('G:/大量データ解析/データ格納庫/BLUE-JAY/m0_BLUE-JAY_by_FFT.csv',"w")
f2.write("year,month,day,time,kn,HDG,H1,theta1,lamda1,w1,H2,theta2,lambda2,w2,error,HbyJWA,TbyJWA,keido,ido,P,R,H,S,")
f2.write('\n')

"""実行部分"""
def Searching(date_number,listname):
    judge=date_number in listname
    if date_number in listname:
        return listname.index(date_number),judge
    return 0,judge
    
def HMS(Date):
    row, judge=Searching(Date,HMS_list)
    if judge==True:
        HMSdata=pd.read_csv(glob.glob(HMS_path)[row],engine='python',usecols=[5,6],header=0).dropna()
        HMSdata=HMSdata.values
        judge=False
        if len(HMSdata) >= fs_HMS*60*50:#当該csvでデータ量を確認
            judge=True
            return judge , np.array(HMSdata[:,0])/2+np.array(HMSdata[:,1])/2
    return judge,0
    
def SIMS(Date):
    row,judge=Searching(Date,SIMS_list)
    if judge==True:
        SIMSdata=pd.read_csv(glob.glob(SIMS_path)[row],engine='python',usecols=[1,2,8],header=0).dropna()
        SIMSdata=SIMSdata.values
        judge=False
        if len(SIMSdata) >= fs_SIMS*60*50:#当該csvでデータ量を確認
            judge=True
            return judge , np.array(SIMSdata[:,0]), np.array(SIMSdata[:,1]), np.array(SIMSdata[:,2])
    return judge, 0,0,0
    
def VDR(Date):
    row , judge=Searching(Date,VDR_list)
    if judge==True:
        VDRdata=pd.read_csv(glob.glob(VDR_path)[row],engine='python',usecols=[7,8],header=0).dropna()#status_water_speed,HDG
        VDRdata=VDRdata.values
        judge=False
        if len(VDRdata) >= fs_VDR*60*50:#当該csvでデータ量を確認
            if (np.max(VDRdata[:,0])-np.min(VDRdata[:,0])) <=2:#1時間で2kn以上の変化がないことを確認
                #角度変化の時系列作成
                direction_t=np.zeros(len(VDRdata))
                for i in range(len(VDRdata)-1):
                    direction_t[i+1] = np.arctan2(np.sin((VDRdata[i+1,1]-VDRdata[i,1])*np.pi/180),np.cos((VDRdata[i+1,1]-VDRdata[i,1])*np.pi/180))*180/np.pi+ direction_t[i]
                #1時間で船首方位が15°以内であることを確認
                if (np.max(direction_t)-np.min(direction_t)) <= 15 :
                    judge=True
                    return judge , np.mean(VDRdata[:,0]), VDRdata[0,1]+np.mean(VDRdata[:,1])
    return judge,0,0

def FFT(time_series,fs):
    freq=np.linspace(0,fs,len(time_series))
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
    f.write(str(Date[0:4])+","+str(Date[4:6])+","+str(Date[6:8])+","+str(Date[8:10])+","+str(speed_mean)+","+str(HDG_mean)+","+
        str(basis[n_basis][0])+","+str(basis[n_basis][1])+","+str(basis[n_basis][2])+","+str(basis[n_basis][3])+","+
        str(basis[n_basis][4])+","+str(basis[n_basis][5])+","+str(basis[n_basis][6])+","+str(basis[n_basis][7])+","+str(basis[n_basis][8])+","+
        str(basis[n_basis][17])+","+str(basis[n_basis][18])+","+str(basis[n_basis][23])+","+str(basis[n_basis][24])+",") 

    
class Calculation:
    def __init__(self,pitch_t,roll_t,heave_t,GMPGMS_t):
        self.pitch_t =np.real(np.fft.ifft(FFT(pitch_t,fs_SIMS)))*len(pitch_t)
        self.roll_t  =np.real(np.fft.ifft(FFT(roll_t,fs_SIMS)))*len(roll_t)
        self.heave_t =self.Heave_calculation(heave_t)
        self.GMPGMS_t=np.real(np.fft.ifft(FFT(GMPGMS_t,fs_HMS)))*len(GMPGMS_t)
        self.encounter_omega=np.arange(0,2.11,0.03)
        self.slide_size=100
        self.frame_size=2**13
        self.freq=0
        
    def Writing(self,x,f):
        X=interpolate.interp1d(self.freq,x,kind='nearest')
        for i in range(len(self.encounter_omega)):
            f.write(str(X(self.encounter_omega)[i])+",")

    def Writing_m0(self,x,f):
        X=interpolate.interp1d(self.freq,x,kind='nearest')
        m0=0
        for i in range(len(self.encounter_omega)-1):
            m0+=(X(self.encounter_omega)[i]+X(self.encounter_omega)[i+1])*0.03*2
        f.write(str(m0)+",")
    
    def Heave_calculation(self,x):#Heave成分の時系列データ作成
        Z_acc=np.real(np.fft.ifft(FFT(x,fs_SIMS)))*len(x)#周波数カット
        Z_vel=Integration(Z_acc,fs_SIMS)#台形積分
        Z_vel=np.real(np.fft.ifft(FFT(Z_vel,fs_SIMS)))*len(Z_vel)#周波数カット
        Z_dis=Integration(Z_vel,fs_SIMS)#台形積分
        Z_dis=np.real(np.fft.ifft(FFT(Z_dis,fs_SIMS)))*len(Z_dis)#周波数カット
        heave_actual=Z_dis-49.55878*np.sin(self.pitch_t*np.pi/180)#ピッチ成分カット
        heave_actual=np.real(np.fft.ifft(FFT(heave_actual,fs_SIMS)))*len(heave_actual)#周波数カット
        """内挿"""
        X=interpolate.interp1d(np.arange(0,len(heave_actual)/fs_HMS,1/fs_HMS),heave_actual,kind='nearest')
        return X(np.arange(0,len(heave_actual)/fs_HMS-10,1/fs_SIMS))
    
    def Welch_ModifiedPeriodogram(self,X,Y):
        freq,csd=signal.csd(X,Y,fs_SIMS,window='hann',nperseg=self.frame_size,noverlap=self.frame_size-self.slide_size,scaling='density')
        self.freq=freq
        return csd/(2*np.pi)

def main():
    print("Start :"+str(len(basis))+"回")
    for n_basis in range(len(basis)):
        Date=str(basis[n_basis][9])+str('%02.0f'%(basis[n_basis][10]))+str('%02.0f'%(basis[n_basis][11]))+str('%02.0f'%(basis[n_basis][12]))
        VDRjudge, speed_mean , HDG_mean = VDR(Date)#VDRデータから対象csvファイル検索
        if VDRjudge == False:
            continue
        
        SIMSjudge, pitch_t, roll_t ,heave_t = SIMS(Date)#SIMSデータから対象csvファイル検索
        if SIMSjudge == False:
            continue
        
        HMSjudge, GMPGMS_t = HMS(Date)#HMSデータから対象csvファイル検索
        if HMSjudge==False:
            continue
        
        DataSet=Calculation(pitch_t,roll_t,heave_t,GMPGMS_t)
        """計算&書き込み開始"""
        print(n_basis)
        Writing_basis(n_basis,Date,speed_mean,HDG_mean,f1)
        Writing_basis(n_basis,Date,speed_mean,HDG_mean,f2)
        """Power Spectrum Density"""
        XY = DataSet.Welch_ModifiedPeriodogram(pitch_t,pitch_t)
        DataSet.Writing(np.abs(XY),f1)
        DataSet.Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(roll_t,roll_t)
        DataSet.Writing(np.abs(XY),f1)
        DataSet.Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(heave_t,heave_t)
        DataSet.Writing(np.abs(XY),f1)
        DataSet.Writing_m0(np.abs(XY),f2)
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,GMPGMS_t)
        DataSet.Writing(np.abs(XY),f1)
        DataSet.Writing_m0(np.abs(XY),f2)
        """Cross Power SPectrum Density"""
        XY = DataSet.Welch_ModifiedPeriodogram(pitch_t,roll_t)
        DataSet.Writing(XY.real,f1)
        DataSet.Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(roll_t,heave_t)
        DataSet.Writing(XY.real,f1)
        DataSet.Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(heave_t,pitch_t)
        DataSet.Writing(XY.real,f1)
        DataSet.Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,roll_t)
        DataSet.Writing(XY.real,f1)
        DataSet.Writing(XY.imag,f1)
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,heave_t)
        DataSet.Writing(XY.real,f1)
        DataSet.Writing(XY.imag,f1)        
        XY = DataSet.Welch_ModifiedPeriodogram(GMPGMS_t,pitch_t)
        DataSet.Writing(XY.real,f1)
        DataSet.Writing(XY.imag,f1)
        
        f1.write('\n')
        f2.write('\n')
    print("End")
    f1.close()
    f2.close()

if __name__ == '__main__':
    main()
