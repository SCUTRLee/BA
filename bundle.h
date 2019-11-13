#ifndef BUNDLE_H_INCLUDED
#define BUNDLE_H_INCLUDED

#include "stdafx.h"
#include "opencv.h"


//�Ƿ���У׼
enum selfParam{
	SELF_CALIB_ON = 0,
	SELF_CALIB_OFF = 1
};

//������У׼
enum selfAlterParam{
	ALTERNATE_ON = 0,
	ALTERNATE_OFF = 1
};

//��������������
enum fixParam{
	//R0��T0��7���̶�
	FIX_7_AUTO = 0, //�Ԅ��ж�
	FIX_7_TX = 1, //tx
	FIX_7_TY = 2, //tx
	FIX_7_TZ = 3, //tx
	//6�̶� R0,T0
	FIX_6 = 4,
	//���̶������в�����������
	FIX_OFF = 5
};

//�Ƿ��������Ľ������Ҫ��
enum fcParam{
	FC_VARIABLE = 0,
	FC_FIX = 1
};

class Bundle
{
private:
	Mat& points; //3D(X,Y,Z)
	Mat& imageX; //ͼ������x
	Mat& imageY; //ͼ������y
	Mat * R; // 3*3
	Mat * T; // 3*1
	Mat * K; // 3*3

	//����ͼƬ�ܵ�����֡����ͼ��ߴ�
	const int pointTotal, frameTotal, width, height;

	//���ռ��еĵ�δͶӰ�����ʱ��ͼ������
	//��ʼ��Ϊ-999.
	float imageNaN;

	//Levenberg-Marquardt����c
	//�ڹ��캯���г�ʼ��Ϊ0.00001F
	float c;

	//��ֹ������
	//�ڹ��캯���г�ʼ��Ϊ0.1
	float epsilon; 

	//�����������
	//�ڹ��캯���г�ʼ��Ϊ5
	//4 : K1,K2,P1,P2
	//5 : K1,K2,K3,P1,P2
	int paramD;

public:
	selfParam selfFlag; //�Ƿ���У׼
	selfAlterParam selfAlterFlag; //������У׼
	fixParam fixFlag; //��������������
	fcParam fcFlag; //�Ƿ��������Ľ������Ҫ��
	float E, newE; //��ͶӰ���E
	int allSize; //��������
	int useAllSize; //Ҫ���µĲ���
	int freedom; //���ɶ�
	int param; //ÿһ��ͼƬ�Ĳ�������
	int paramK; //�ڲ� 3��3 + paraD��
	bool cFlag; //C��־���ڵ�������
	int KPtotal; //����ͼƬ�������������ܺͣ���������λ���깲��509������KPtotal = 509 * 5
	Mat Distort; //����ϵ������
	vector<float> vectorE; //���ڴ洢��ͶӰ���E������


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

		//�����Ͻ�����ת��Ϊ��������
		changeCoordToCenter();

		//fx = fy =��fx + fy) / 2
		fxfyAverage();

		//����KPtotal
		calcKPtotal();
		
		cout << "Bundle::Constructor" << endl;
		cout << "frameTotal = " << frameTotal << endl;
		cout << "pointTotal = " << pointTotal << endl;
		cout << "KPTotal = " << KPtotal << endl << endl;

	}

	///-----BA��ʼ-----
	void start(selfParam _selfFlag=SELF_CALIB_ON, fixParam _fixFlag=FIX_7_AUTO, fcParam _fcFlag=FC_VARIABLE, selfAlterParam _selfAlterFlag=ALTERNATE_ON);

	float calcE();

	///-----����H,g-----
	void calcHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf);

	///-----����һ�η���ʽ�ļ���-----
	void solveHg(Mat& HE,Mat& HF,Mat& HG,Mat& Gp,Mat& Gf,Mat& Dp,Mat& Df);

	///-----����k�η���ʽ�ļ���-----
	void solveHg_k(Mat& HG,Mat& Gf,Mat& Dk);

	///-----���½�-----
	void update(Mat& Dp,Mat& Df);

	///-----���½�-----
	void update_6(Mat& Dp,Mat& Df);

	///-----���½�-----
	void update_k(Mat& Dk);

	///-----���ؽ������-----
	void unupdate(Mat& Dp,Mat& Df);

	///-----���ؽ������-----
	void unupdate_6(Mat& Dp,Mat& Df);

	///-----���ؽ������-----
	void unupdate_k(Mat& Dk);

	///-----�ƶ��ض����к���-----
	void moveLine(Mat& A,int i_=-1,int i=-1,int j_=-1,int j=-1);

	///-----�����Ͻ�����ת��Ϊ��������-----
	void changeCoordToCenter();

	///-----����������ת��Ϊ��������-----
	void changeCoordToLeftUpper();

	///-----ͳһ���࣬��fx��fy��f-----
	void fxfyAverage();

	///-----����KPtotal-----
	void calcKPtotal();

	///-----imageNaN���-----
	//imageNaN��ʼֵ-999.
	void set_imageNaN(float _imageNaN)
	{
		cout << "Set imageNaN " << imageNaN;
		imageNaN = _imageNaN;
		cout << " to " << imageNaN << endl << endl;
	}

	///-----c���-----
	//c��ʼֵ0.00001F
	void set_c(float _c)
	{
		cout << "Set c " << c;
		c = _c;
		cout << " to " << c << endl << endl;
	}

	///-----epsilon���-----
	//epsilon��ʼֵ0.1
	void set_epsilon(float _epsilon)
	{
		cout << "Set epsilon " << epsilon;
		epsilon = _epsilon;
		cout << " to " << epsilon << endl << endl;
	}

	///-----paramD���-----
	//paramD��ʼֵ5
	void set_paramD(int _paramD)
	{
		cout << "Set paramD " << paramD;
		paramD = _paramD;
		cout << " to " << paramD << endl << endl;

		//���������ʼ��
		Distort = (Mat_<float>(paramD,frameTotal));
		Distort = Scalar::all(0);
	}

	///-----Distort����-----
	void set_Distort(Mat& _Distort)
	{
		cout << "Set Distort " << endl;

		Distort = _Distort;
	}

	///-----�������΢��-----
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
	
	///-----һ��΢�ּ���-----
	inline float dE(float dp, float dq, float dr, float p, float q, float r, float x0, float y0, float dx=0, float dy=0, float ddx=0, float ddy=0)
	{
		return 2*((p/r+dx-x0)*((r*dp-p*dr)/(r*r)+ddx)+(q/r+dy-y0)*((r*dq-q*dr)/(r*r)+ddy));
	}
	
	///-----����΢�ּ���-----
	inline float ddE(float dpk, float dqk, float drk, float dpl, float dql, float drl, float p, float q, float r, float ddxk=0, float ddyk=0, float ddxl=0, float ddyl=0)
	{
		return 2*(((r*dpk-p*drk)/(r*r)+ddxk)*((r*dpl-p*drl)/(r*r)+ddxl)+((r*dqk-q*drk)/(r*r)+ddyk)*((r*dql-q*drl)/(r*r)+ddyl));
	}

};

#endif