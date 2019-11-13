#include "stdafx.h"
#include "opencv.h" 
#include "file.h" 
#include "bundle.h" 

//全局变量声明

selfParam selfFlag = SELF_CALIB_ON;
selfAlterParam selfAlterFlag = ALTERNATE_OFF;
fixParam fixFlag = FIX_7_AUTO;
fcParam fcFlag = FC_FIX;
float c = 0.0001F, epsilon = 0.01F;
const int FRAMETOTAL = 5; //帧数
const int BNDL = FRAMETOTAL - 1; //捆束宽度（帧数-1
const int WID = 128, HEI = 128; //图像大小
const char filedir[] = "data"; //放置数据的目录名称
const char currentdir[] = "sample"; // 要使用的目录名称-它最终将以./filedir/currentdir/hoge.txt的形式出现

const char corrfilename[] = "allbundlepoints.csv"; //特征点输入数据
const char rtfilename[] = "rt.xml"; //平移/旋转矢量输入数据
const char camerafilename[] = "camera.xml"; //相机内部矩阵输入数据
const char corroutfilename[] = "result_xyz"; //特征点输出数据
const char rtoutfilename[] = "result_rt.xml"; //平移/旋转矢量输出
const char Routfilename[] = "result_R.xml"; //平移/旋转矢量输出
const char Toutfilename[] = "result_T.xml"; //平移/旋转矢量输出
const char cameraoutfilename[] = "result_camera.xml"; //摄像机内部矩阵输出
const char distortoutfilename[] = "result_distort.txt"; //摄像机畸变系数输出
int KPtotal; //特征点的实际数量“ KP =关键点=特征点”
int KP2Dtotal; //每帧中的特征点总数
float fx, fy, cx, cy; //每帧特征点的数量总浮点数fx，fy，cx，cy
//double PI = 3.14159265358979;
//int fps = 15; 

////快速排序比较功能
int comp( const void *c1, const void *c2 );


//////////////////
// main函数
//
int _tmain(int argc, _TCHAR* argv[])
{
	cout << "-------------------------" << endl;
	cout << "    Bundle Adjustment    " << endl;
	cout << "-------------------------" << endl;
	clock_t start_time_total,end_time_total;
	start_time_total = clock(); //执行时间测量


	///---特征点输入数据读取---
	char filename[255];
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,corrfilename);

	// cout <<“特征点输入数据” <<文件名<<“ Read” << endl; 
	//获取文件行数（=特征点数）并分配给KPtotal 
	KPtotal = readfileLine(filename); // readfileLine函数-请参见file.cpp
	// cout <<“完成的特征点数量=” << KPtotal <<“件” << endl << endl; 
	//放置特征点的垫子放置图像坐标的垫子

	Mat points = Mat::zeros(KPtotal,3,CV_32F);
	Mat imageX = Mat::zeros(KPtotal,FRAMETOTAL,CV_32F); // 每张图的图像坐标X
	Mat imageY = Mat::zeros(KPtotal,FRAMETOTAL,CV_32F); // 每张图的图像坐标Y

	//读取文件并分配给Mat
	readfilePoints(filename, points, imageX, imageY);
	
	//结构CorrPoint实例化结构声明在file.h中
	CorrPoint * corrpoint;
	corrpoint = (CorrPoint *)malloc(sizeof(CorrPoint) * KPtotal);
	initCorrPoint(corrpoint); //对应点初始化-请参阅file.cpp
    //读取文件并分配给corrpoint
	readfileCorr(filename, corrpoint); // readfileCorr函数-请参见file.cpp
	


	// ---摄像机平移/旋转矢量读取---
	// 当使用ICPnormalMultiData输出rt.xml时，请注意RT矩阵！它与正常的4 * 4RT矩阵形式不同。
	// 第四行是[tx，ty，tz，1]，而不是[0,0,0,1]。
	Mat rt[FRAMETOTAL]; //创建矩阵以放置RT4 * 4矩阵

	for(int i=0;i<FRAMETOTAL;i++) 
		rt[i] = Mat_<float>(4,4); // 4 * 4内存分配

	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,rtfilename);
	readfileRT(filename, rt);


	/// ---摄像机内部矩阵（fx，fy，cx，cy）读取---
	////通过OpenCV摄像机校准获得的（Cx，cy）位于图像的左上方中心，
	////此程序 然后，（cx，cy）表示距图像中心的偏差。
	////稍微向前，转换为从中心偏移
	////在这里，我们只读取数据
	Mat cameraA(3,3,CV_32F); //相机内部矩阵3 * 3矩阵 [fx, 0, cx; 0, fy, cy; 0, 0, 1]
	//读取camera.xml文件并分配给Mat
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,camerafilename);
	//cout << "相机内部矩阵数据" << filename << " 已读" << endl;
	readfileCamera(filename, &cameraA); //readfileCamera函数-请参见file.cpp
	//cout << "完成" << endl << endl;


	///---Bundle数据格式传递---
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
	//cout << "   开始BA   " << endl;
	//cout << "------------------------------" << endl;
	clock_t start_time_bundle,end_time_bundle;
	start_time_bundle = clock();


	Bundle bundle(points, imageX, imageY, R, T, K, WID, HEI);

	///---条件设定---
	//bundle.set_imageNaN(-1.); ///未投影在相机上时的图像坐标

	bundle.set_c(c);     //Levenberg-Marquardt方法cc
	bundle.set_epsilon(epsilon);    //终止条件ε

	///---BA-start---
	bundle.start(selfFlag, fixFlag, fcFlag, selfAlterFlag);
	

	end_time_bundle = clock(); //执行时间结束
	cout << "总时长 = " << (float)(end_time_bundle - start_time_bundle)/CLOCKS_PER_SEC << "秒" << endl << endl;



	
	///---点集(x,y,z)输出---
	///(x,y,z)的CSV文件
	sprintf_s(filename, _countof(filename), "./%s/%s/%s.csv", filedir, currentdir,corroutfilename);
	writefilePoints(filename, points, imageX.rows); 
	///PointCloudLibrary点集格式.pcd
	//根据z的值，颜色指定为 近-远 绿-红
	//求出z的最大值和最小值
	double minVal, maxVal;
	Point minLoc, maxLoc;
	minMaxLoc(points(Range(0,points.rows),Range(2,3)), &minVal, &maxVal, &minLoc, &maxLoc);
	//写入文件
	sprintf_s(filename, _countof(filename), "./%s/%s/%s.pcd", filedir, currentdir,corroutfilename);
	writefilePoints(filename, points, imageX.rows, imageX.cols, (float)minVal, (float)maxVal, R, T); //writefile娭悢 - file.cpp嶲徠
	

	///---夞揮R偺弌椡---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,Routfilename);
	writefileMat(filename, R); //writefileMat娭悢 - file.cpp嶲徠
	

	///---暲恑T偺弌椡---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,Toutfilename);
	writefileMat(filename, T); //writefileMat娭悢 - file.cpp嶲徠

	
	///---僇儊儔撪晹峴楍(fx,fy,cx,cy)偺弌椡---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,cameraoutfilename);
	writefileMat(filename, K); //writefileMat娭悢 - file.cpp嶲徠

	
	///---僇儊儔榗傒學悢(K1,K2,P1,P2,K3)偺弌椡---
	sprintf_s(filename, _countof(filename), "./%s/%s/%s", filedir, currentdir,distortoutfilename);
	writefileDistort(filename, bundle.Distort); //writefileMat娭悢 - file.cpp嶲徠

	
	///---嵞搳塭岆嵎vectorE彂偒崬傒---
	FILE * fp;
	sprintf_s(filename, _countof(filename), "%s/%s/result_E.csv",filedir,currentdir);
	errno_t error;
	error = fopen_s(&fp, filename, "w");
	if(error != 0){
		cout << "僼傽僀儖偑奐偗傑偣傫 " << filename << endl;
		exit(1);
	}
	for(int i=0;i<bundle.vectorE.size();i++){
		fprintf(fp,"%f\n",bundle.vectorE[i]);
	}
	fclose(fp);
	
	
	///---僾儘僌儔儉廔椆張棟---
	cout << endl;
	cout << "--------------" << endl;
	cout << "    Finish    " << endl;
	cout << "--------------" << endl;
	end_time_total = clock(); //幚峴帪娫寁應廔椆
	cout << "僾儘僌儔儉幚峴帪娫 = " << (float)(end_time_total - start_time_total)/CLOCKS_PER_SEC << "昩" << endl << endl;


	return 0;
}



///////////////////////////////
// 僋僀僢僋僜乕僩梡斾妑娭悢
//
int comp( const void *c1, const void *c2 )
{
	CorrPoint point1 = *(CorrPoint *)c1;
	CorrPoint point2 = *(CorrPoint *)c2;

	float tmp1 = point1.z;   /* z 傪婎弨偲偡傞 */
	float tmp2 = point2.z;

	if(tmp1 == tmp2)       return 0;
	else if(tmp1 > tmp2)   return 1;
	else					 return -1;
}