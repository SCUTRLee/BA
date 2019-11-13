#ifndef BUNDLE_H_INCLUDED
#define BUNDLE_H_INCLUDED

#include "stdafx.h"
#include "opencv.h"


//是否自校准
enum selfParam{
	SELF_CALIB_ON = 0,
	SELF_CALIB_OFF = 1
};

//交替自校准
enum selfAlterParam{
	ALTERNATE_ON = 0,
	ALTERNATE_OFF = 1
};

//修正参数的设置
enum fixParam{
	//R0，T0和7个固定
	FIX_7_AUTO = 0, //自動判定
	FIX_7_TX = 1, //tx
	FIX_7_TY = 2, //tx
	FIX_7_TZ = 3, //tx
	//6固定 R0,T0
	FIX_6 = 4,
	//不固定，所有参数可以修正
	FIX_OFF = 5
};

//是否更新相机的焦距和主要点
enum fcParam{
	FC_VARIABLE = 0,
	FC_FIX = 1
};

class Bundle
{
private:
	Mat& points; //3D(X,Y,Z)
	Mat& imageX; //图像坐标x
	Mat& imageY; //图像坐标y
	Mat * R; // 3*3
	Mat * T; // 3*1
	Mat * K; // 3*3

	//所有图片总点数，帧数，图像尺寸
	const int pointTotal, frameTotal, width, height;

	//当空间中的点未投影到相机时的图像坐标
	//初始化为-999.
	float imageNaN;

	//Levenberg-Marquardt方法c
	//在构造函数中初始化为0.00001F
	float c;

	//终止条件ε
	//在构造函数中初始化为0.1
	float epsilon; 

	//畸变参数数量
	//在构造函数中初始化为5
	//4 : K1,K2,P1,P2
	//5 : K1,K2,K3,P1,P2
	int paramD;

public:
	selfParam selfFlag; //是否自校准
	selfAlterParam selfAlterFlag; //交替自校准
	fixParam fixFlag; //修正参数的设置
	fcParam fcFlag; //是否更新相机的焦距和主要点
	float E, newE; //重投影误差E
	int allSize; //参数总数
	int useAllSize; //要更新的参数
	int freedom; //自由度
	int param; //每一张图片的参数数量
	int paramK; //内参 3或（3 + paraD）
	bool cFlag; //C标志用于迭代计算
	int KPtotal; //所有图片中特征点数的总和，这里面三位坐标共有509个，即KPtotal = 509 * 5
	Mat Distort; //畸变系数矩阵
	vector<float> vectorE; //用于存储重投影误差E的向量


	Bundle(Mat& _points, Mat& _imageX, Mat& _imageY, 
		Mat * _R, Mat * _T, Mat * _K, int _width, int _height) : 
		points(_points), imageX(_imageX), imageY(_imageY), width(_width), height(_height), 
		imageNaN(-999.), epsilon(0.1F), c(0.0001F), 
		pointTotal(imageX.rows), frameTotal(imageX.cols), paramD(5)
	{
		R = _R;
		T = _T;
		K = _K;

		Distort = (Mat_<float>(paramD,frameTotal));
		Distort = Scalar::all(0);

		//从左上角坐标转换为中心坐标
		changeCoordToCenter();

		//fx = fy =（fx + fy) / 2
		fxfyAverage();

		//计算KPtotal
		calcKPtotal();
		
		cout << "Bundle::Constructor" << endl;
		cout << "frameTotal = " << frameTotal << endl;
		cout << "pointTotal = " << pointTotal << endl;
		cout << "KPTotal = " << KPtotal << endl << endl;

	}

	///-----BA开始-----
	void start(selfParam _selfFlag=SELF_CALIB_ON, fixParam _fixFlag=FIX_7_AUTO, fcParam _fcFlag=FC_VARIABLE, selfAlterParam _selfAlterFlag=ALTERNATE_ON);

	float calcE();

	///-----计算H,g-----
	void calcHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf);

	///-----联立一次方程式的计算-----
	void solveHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf,Mat& Dp,Mat& Df);

	///-----联立k次方程式的计算-----
	void solveHg_k(Mat& HG,Mat& Gf,Mat& Dk);

	///-----更新解-----
	void update(Mat& Dp,Mat& Df);

	///-----更新解-----
	void update_6(Mat& Dp,Mat& Df);

	///-----更新解-----
	void update_k(Mat& Dk);

	///-----返回解决方案-----
	void unupdate(Mat& Dp,Mat& Df);

	///-----返回解决方案-----
	void unupdate_6(Mat& Dp,Mat& Df);

	///-----返回解决方案-----
	void unupdate_k(Mat& Dk);

	///-----移动特定的行和列-----
	void moveLine(Mat& A,int i_=-1,int i=-1,int j_=-1,int j=-1);

	///-----从左上角坐标转换为中心坐标-----
	void changeCoordToCenter();

	///-----将中心坐标转换为左上坐标-----
	void changeCoordToLeftUpper();

	///-----统一焦距，从fx，fy到f-----
	void fxfyAverage();

	///-----计算KPtotal-----
	void calcKPtotal();

	///-----imageNaN变更-----
	//imageNaN初始值-999.
	void set_imageNaN(float _imageNaN)
	{
		cout << "Set imageNaN " << imageNaN;
		imageNaN = _imageNaN;
		cout << " to " << imageNaN << endl << endl;
	}

	///-----c变更-----
	//c初始值0.00001F
	void set_c(float _c)
	{
		cout << "Set c " << c;
		c = _c;
		cout << " to " << c << endl << endl;
	}

	///-----epsilon变更-----
	//epsilon初始值0.1
	void set_epsilon(float _epsilon)
	{
		cout << "Set epsilon " << epsilon;
		epsilon = _epsilon;
		cout << " to " << epsilon << endl << endl;
	}

	///-----paramD变更-----
	//paramD初始值5
	void set_paramD(int _paramD)
	{
		cout << "Set paramD " << paramD;
		paramD = _paramD;
		cout << " to " << paramD << endl << endl;

		//畸变参数初始化
		Distort = (Mat_<float>(paramD,frameTotal));
		Distort = Scalar::all(0);
	}

	///-----Distort输入-----
	void set_Distort(Mat& _Distort)
	{
		cout << "Set Distort " << endl;

		Distort = _Distort;
	}

	///-----畸变参数微分-----
	inline float calcddx(float dp, float dq, float dr, float p, float q, float r, float cx, float cy, float K1, float K2, float K3, float P1, float P2, float xc, float yc, float r2, float r4, float r6, int cxFlag, int cyFlag)
	{
		float x_ = (r*dp-p*dr)/(r*r);
		float y_ = (r*dq-q*dr)/(r*r);
		float r2_; 
		if(cxFlag==0 && cyFlag==0){
			r2_ = 2*xc*x_ + 2*yc*y_;
		}else if(cxFlag!=0){
			r2_ = 2*xc*(x_-1) + 2*yc*y_;
		}else if(cyFlag!=0){
			r2_ = 2*xc*x_ + 2*yc*(y_-1);
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}
		float r4_ = 2*r2*r2_;
		float r6_ = 3*r4*r2_;

		float ddx;
		if(cxFlag==0 && cyFlag==0){
			ddx = (K1*r2+K2*r4+K3*r6)*x_ + (K1*r2_+K2*r4_+K3*r6_)*xc + P1*(r2_+4*xc*x_) + 2*P2*(x_*yc+xc*y_);
		}else if(cxFlag!=0){
			ddx = (K1*r2+K2*r4+K3*r6)*(x_-1) + (K1*r2_+K2*r4_+K3*r6_)*xc + P1*(r2_+4*xc*(x_-1)) + 2*P2*((x_-1)*yc+xc*y_);
		}else if(cyFlag!=0){
			ddx = (K1*r2+K2*r4+K3*r6)*x_ + (K1*r2_+K2*r4_+K3*r6_)*xc + P1*(r2_+4*xc*x_) + 2*P2*(x_*yc+xc*(y_-1));
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}

		return ddx;
	}

	inline float calcddy(float dp, float dq, float dr, float p, float q, float r, float cx, float cy, float K1, float K2, float K3, float P1, float P2, float xc, float yc, float r2, float r4, float r6, int cxFlag, int cyFlag)
	{
		float x_ = (r*dp-p*dr)/(r*r);
		float y_ = (r*dq-q*dr)/(r*r);
		float r2_; 
		if(cxFlag==0 && cyFlag==0){
			r2_ = 2*xc*x_ + 2*yc*y_;
		}else if(cxFlag!=0){
			r2_ = 2*xc*(x_-1) + 2*yc*y_;
		}else if(cyFlag!=0){
			r2_ = 2*xc*x_ + 2*yc*(y_-1);
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}
		float r4_ = 2*r2*r2_;
		float r6_ = 3*r4*r2_;

		float ddy;
		if(cxFlag==0 && cyFlag==0){
			ddy = (K1*r2+K2*r4+K3*r6)*y_ + (K1*r2_+K2*r4_+K3*r6_)*yc + 2*P1*(x_*yc+xc*y_) + P2*(r2_+4*yc*y_);
		}else if(cxFlag!=0){
			ddy = (K1*r2+K2*r4+K3*r6)*y_ + (K1*r2_+K2*r4_+K3*r6_)*yc + 2*P1*((x_-1)*yc+xc*y_) + P2*(r2_+4*yc*y_);
		}else if(cyFlag!=0){
			ddy = (K1*r2+K2*r4+K3*r6)*(y_-1) + (K1*r2_+K2*r4_+K3*r6_)*yc + 2*P1*(x_*yc+xc*(y_-1)) + P2*(r2_+4*yc*(y_-1));
		}else{
			cout << "Error at dE" << endl;
			exit(1);
		}

		return ddy;
	}
	
	///-----一次微分计算-----
	inline float dE(float dp, float dq, float dr, float p, float q, float r, float x0, float y0, float dx=0, float dy=0, float ddx=0, float ddy=0)
	{
		return 2*((p/r+dx-x0)*((r*dp-p*dr)/(r*r)+ddx)+(q/r+dy-y0)*((r*dq-q*dr)/(r*r)+ddy));
	}
	
	///-----二次微分计算-----
	inline float ddE(float dpk, float dqk, float drk, float dpl, float dql, float drl, float p, float q, float r, float ddxk=0, float ddyk=0, float ddxl=0, float ddyl=0)
	{
		return 2*(((r*dpk-p*drk)/(r*r)+ddxk)*((r*dpl-p*drl)/(r*r)+ddxl)+((r*dqk-q*drk)/(r*r)+ddyk)*((r*dql-q*drl)/(r*r)+ddyl));
	}

};

#endif