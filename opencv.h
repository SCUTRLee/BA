#ifndef OPENCV_INCLUDED
#define OPENCV_INCLUDED

////项目属性⇒C \\ C ++常规
//// 使用opencv2添加文件夹“ C：\ OpenCV \ include”等。

#include "opencv2\opencv.hpp"

using namespace cv;
#ifdef _DEBUG
    ////Debug
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_core249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_imgproc249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_highgui249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_objdetect249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_contrib249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_features2d249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_flann249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_gpu249d.lib")
    //#pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_haartraining_engined.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_legacy249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_ts249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_video249d.lib")
	////捛壛偟偨傕偺
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_ml249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_nonfree249d.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_calib3d249d.lib")

#else
    ////Release
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_core249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_imgproc249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_highgui249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_objdetect249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_contrib249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_features2d249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_flann249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_gpu249.lib")
    //#pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_haartraining_engined.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_legacy249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_ts249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_video249.lib")
	////捛壛偟偨傕偺
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_ml249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_nonfree249.lib")
    #pragma comment(lib,"E:\\opencv\\2.4.9\\opencv\\build\\x86\\vc10\\lib\\opencv_calib3d249.lib")
#endif

#endif