#include "stdafx.h"
#include "opencv.h" 
#include "file.h" 
#include "bundle.h" 

//ȫ�ֱ�������

selfParam selfFlag = SELF_CALIB_ON;
selfAlterParam selfAlterFlag = ALTERNATE_OFF;
fixParam fixFlag = FIX_7_AUTO;
fcParam fcFlag = FC_FIX;
float c = 0.0001F, epsilon = 0.01F;
const int FRAMETOTAL = 5; //֡��
const int BNDL = FRAMETOTAL - 1; //������ȣ�֡��-1
const int WID = 128, HEI = 128; //ͼ���С
const char filedir[] = "data"; //�������ݵ�Ŀ¼����
const char currentdir[] = "sample"; // Ҫʹ�õ�Ŀ¼����-�����ս���./filedir/currentdir/hoge.txt����ʽ����

const char corrfilename[] = "allbundlepoints.csv"; //��������������
const char rtfilename[] = "rt.xml"; //ƽ��/��תʸ����������
const char camerafilename[] = "camera.xml"; //����ڲ�������������
const char corroutfilename[] = "result_xyz"; //�������������
const char rtoutfilename[] = "result_rt.xml"; //ƽ��/��תʸ�����
const char Routfilename[] = "result_R.xml"; //ƽ��/��תʸ�����
const char Toutfilename[] = "result_T.xml"; //ƽ��/��תʸ�����
const char cameraoutfilename[] = "result_camera.xml"; //������ڲ��������
const char distortoutfilename[] = "result_distort.txt"; //���������ϵ�����
int KPtotal; //�������ʵ�������� KP =�ؼ���=�����㡱
int KP2Dtotal; //ÿ֡�е�����������
float fx, fy, cx, cy; //ÿ֡������������ܸ�����fx��fy��cx��cy
//double PI = 3.14159265358979;
//int fps = 15; 

////��������ȽϹ���
int comp( const void *c1, const void *c2 );


//////////////////
// main����
//
int _tmain(int argc, _TCHAR* argv[])
{
	cout << "-------------------------" << endl;
	cout << "    Bundle Adjustment    " << endl;
	cout << "-------------------------" << endl;
	clock_t start_time_total,end_time_total;
	start_time_total = clock(); //ִ��ʱ�����


	///---�������������ݶ�ȡ---
	char filename[255];
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,corrfilename);

	// cout <<���������������ݡ� <<�ļ���<<�� Read�� << endl; 
	//��ȡ�ļ�������=�����������������KPtotal 
	KPtotal = readfileLine(filename); // readfileLine����-��μ�file.cpp
	// cout <<����ɵ�����������=�� << KPtotal <<������ << endl << endl; 
	//����������ĵ��ӷ���ͼ������ĵ���

	Mat points = Mat::zeros(KPtotal,3,CV_32F);
	Mat imageX = Mat::zeros(KPtotal,FRAMETOTAL,CV_32F); // ÿ��ͼ��ͼ������X
	Mat imageY = Mat::zeros(KPtotal,FRAMETOTAL,CV_32F); // ÿ��ͼ��ͼ������Y

	//��ȡ�ļ��������Mat
	readfilePoints(filename, points, imageX, imageY);
	
	//�ṹCorrPointʵ�����ṹ������file.h��
	CorrPoint * corrpoint;
	corrpoint = (CorrPoint *)malloc(sizeof(CorrPoint) * KPtotal);
	initCorrPoint(corrpoint); //��Ӧ���ʼ��-�����file.cpp
    //��ȡ�ļ��������corrpoint
	readfileCorr(filename, corrpoint); // readfileCorr����-��μ�file.cpp
	


	// ---�����ƽ��/��תʸ����ȡ---
	// ��ʹ��ICPnormalMultiData���rt.xmlʱ����ע��RT��������������4 * 4RT������ʽ��ͬ��
	// ��������[tx��ty��tz��1]��������[0,0,0,1]��
	Mat rt[FRAMETOTAL]; //���������Է���RT4 * 4����

	for(int i=0;i<FRAMETOTAL;i++) 
		rt[i] = Mat_<float>(4,4); // 4 * 4�ڴ����

	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,rtfilename);
	readfileRT(filename, rt);


	/// ---������ڲ�����fx��fy��cx��cy����ȡ---
	////ͨ��OpenCV�����У׼��õģ�Cx��cy��λ��ͼ������Ϸ����ģ�
	////�˳��� Ȼ�󣬣�cx��cy����ʾ��ͼ�����ĵ�ƫ�
	////��΢��ǰ��ת��Ϊ������ƫ��
	////���������ֻ��ȡ����
	Mat cameraA(3,3,CV_32F); //����ڲ�����3 * 3���� [fx, 0, cx; 0, fy, cy; 0, 0, 1]
	//��ȡcamera.xml�ļ��������Mat
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,camerafilename);
	//cout << "����ڲ���������" << filename << " �Ѷ�" << endl;
	readfileCamera(filename, &cameraA); //readfileCamera����-��μ�file.cpp
	//cout << "���" << endl << endl;


	///---Bundle���ݸ�ʽ����---
	Mat R[FRAMETOTAL]; //3*3
	Mat T[FRAMETOTAL]; //3*1
	Mat K[FRAMETOTAL]; //3*3
	for(int i=0;i<FRAMETOTAL;i++){
		R[i] = Mat_<float>(3,3); 
		R[i] = rt[i](Range(0,3),Range(0,3));
		T[i] = Mat_<float>(3,1); 
		T[i] = rt[i](Range(3,4),Range(0,3)).t();
		K[i] = Mat_<float>(3,3); 
		cameraA.copyTo(K[i]); 
	}


	///---BA---
	//cout << endl;
	//cout << "------------------------------" << endl;
	//cout << "   ��ʼBA   " << endl;
	//cout << "------------------------------" << endl;
	clock_t start_time_bundle,end_time_bundle;
	start_time_bundle = clock();


	Bundle bundle(points, imageX, imageY, R, T, K, WID, HEI);

	///---�����趨---
	//bundle.set_imageNaN(-1.); ///δͶӰ�������ʱ��ͼ������

	bundle.set_c(c);     //Levenberg-Marquardt����cc
	bundle.set_epsilon(epsilon);    //��ֹ������

	///---BA-start---
	bundle.start(selfFlag, fixFlag, fcFlag, selfAlterFlag);
	

	end_time_bundle = clock(); //ִ��ʱ�����
	cout << "��ʱ�� = " << (float)(end_time_bundle - start_time_bundle)/CLOCKS_PER_SEC << "��" << endl << endl;



	
	///---�㼯(x,y,z)���---
	///(x,y,z)��CSV�ļ�
	sprintf_s(filename, _countof(filename), "./%s/%s/%s.csv", filedir, currentdir,corroutfilename);
	writefilePoints(filename, points, imageX.rows); 
	///PointCloudLibrary�㼯��ʽ.pcd
	//����z��ֵ����ɫָ��Ϊ ��-Զ ��-��
	//���z�����ֵ����Сֵ
	double minVal, maxVal;
	Point minLoc, maxLoc;
	minMaxLoc(points(Range(0,points.rows),Range(2,3)), &minVal, &maxVal, &minLoc, &maxLoc);
	//д���ļ�
	sprintf_s(filename, _countof(filename), "./%s/%s/%s.pcd", filedir, currentdir,corroutfilename);
	writefilePoints(filename, points, imageX.rows, imageX.cols, (float)minVal, (float)maxVal, R, T); //writefile�֐� - file.cpp�Q��
	

	///---��]R�̏o��---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,Routfilename);
	writefileMat(filename, R); //writefileMat�֐� - file.cpp�Q��
	

	///---���iT�̏o��---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,Toutfilename);
	writefileMat(filename, T); //writefileMat�֐� - file.cpp�Q��

	
	///---�J���������s��(fx,fy,cx,cy)�̏o��---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,cameraoutfilename);
	writefileMat(filename, K); //writefileMat�֐� - file.cpp�Q��

	
	///---�J�����c�݌W��(K1,K2,P1,P2,K3)�̏o��---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,distortoutfilename);
	writefileDistort(filename, bundle.Distort); //writefileMat�֐� - file.cpp�Q��

	
	///---�ē��e�덷vectorE��������---
	FILE * fp;
	sprintf_s(filename, _countof(filename), "%s/%s/result_E.csv",filedir,currentdir);
	errno_t error;
	error = fopen_s(&fp, filename, "w");
	if(error != 0){
		cout << "�t�@�C�����J���܂��� " << filename << endl;
		exit(1);
	}
	for(int i=0;i<bundle.vectorE.size();i++){
		fprintf(fp,"%f\n",bundle.vectorE[i]);
	}
	fclose(fp);
	
	
	///---�v���O�����I������---
	cout << endl;
	cout << "--------------" << endl;
	cout << "    Finish    " << endl;
	cout << "--------------" << endl;
	end_time_total = clock(); //���s���Ԍv���I��
	cout << "�v���O�������s���� = " << (float)(end_time_total - start_time_total)/CLOCKS_PER_SEC << "�b" << endl << endl;


	return 0;
}



///////////////////////////////
// �N�C�b�N�\�[�g�p��r�֐�
//
int comp( const void *c1, const void *c2 )
{
	CorrPoint point1 = *(CorrPoint *)c1;
	CorrPoint point2 = *(CorrPoint *)c2;

	float tmp1 = point1.z;   /* z ����Ƃ��� */
	float tmp2 = point2.z;

	if(tmp1 == tmp2)       return 0;
	else if(tmp1 > tmp2)   return 1;
	else					 return -1;
}