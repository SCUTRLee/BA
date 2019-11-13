#ifndef FILE_H_INCLUDED
#define FILE_H_INCLUDED

#include "stdafx.h"
#include "opencv.h" //OpenCV���C�u����

extern const int FRAMETOTAL,WID,HEI;
extern int KPtotal;

///---���̓f�[�^�����邽�߂̍\����CorrPoint---
//CorrPoint�Ƃ́uCorresponding Point = �Ή��_�v�̂��Ƃł��������_�Ɠ������Ǝv���Ă��������B�����t���[���ԂőΉ����Ă���_�Ȃ̂őΉ��_�ƌĂ�ł��܂�
typedef struct{
	float x;
	float y;
	float z;
	float * x0;
	float * y0;
}CorrPoint;

///---write�֐��Ŏg�p����񋓌^---
enum WriteType{
	XYZ,
	R,
	T,
	FXFYCXCY
};

///---�֐��v���g�^�C�v�錾---
void readfilePoints(char * filename, Mat& points, Mat& imageX, Mat& imageY);
void initCorrPoint(CorrPoint * corrpoint);
void freeCorrPoint(CorrPoint * corrpoint);
void readfileCorr(char * filename, CorrPoint * corrpoint);
int readfileLine(char * filename);
void readfileRT(char * filename, Mat * rt);
void readfileCamera(char * filename, Mat * cameraA);
void writefile(char * filename, float * paraA, WriteType type, float min=-1., float max=-1., Mat * rt = NULL, CorrPoint * corrpoint = NULL);
void writefileMat(char * filename, Mat * rt);
void writefileDistort(char * filename, Mat Distort);
void writefilePoints(char * filename, Mat& points, int pointTotal);
void writefilePoints(char * filename, Mat& points, int pointTotal, int frameTotal, float min, float max, Mat * R = NULL, Mat * T = NULL);

#endif