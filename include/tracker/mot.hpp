/*
 * mot.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: luzdora
 */

#ifndef MOT_HPP_
#define MOT_HPP_

#include "clustersFuncs.hpp"
// LLAMADO A BOOST #include "jmath/jblas.hpp"
#include "klt.h"
#include "klt.hpp"
#include "regions.hpp"
#include "imageConversion.hpp"
#include "observeModel.hpp"
#include "predictModel.hpp"
#include "kalmanFilter.hpp"

namespace tracker {

using namespace boost_incl;
using namespace filter;
enum rect {x1p=0,y1p,x2p,y2p};

class TrackedObject
{
public:

	mat fl;
	mat conv;
	mat contour;

	int curFrame; //! the current frame number in the feature table (frame modulo 96)
	int flagRegion; // 1 :activa el uso de contornos, 0: los puntos y el rectangulo
	int status; // 1: Initialisation 0: not initialisation
	double dt;
	// state
	int scope, prevScope; //! longer distance from one feature to
	//the barycenter

	int boundingBox[4]; //! upperleft and downright corners of the
	//bounding box of the object (x1,y1,x2,y2)

	double pos[2]; //! position of the object in pixels in the image (x,y)
	double vel[2];
	double preVel[2]; //!  previous velocity of the  object in pixels in the image (x,y)
	int maxDisplacement; //! the maximum displacement in pixels from
	//one frame to the next, determining the
	//size of the processing area

	int maxAccel; //! the maximum acceleration of the object in
	//pixels by frames^2, determining the error of the
	//kalman model (speed is constant)

	int measureConfidence; //! the confidence in the measure
	//(barycenter of features representing the object)

	int securityMargin;

	int initSpeedConfidence; //! the confidence in the initial speed,
	// that's to say the maximum value of initial speed
	// (because we assume it to be 0 for the first frame)

	int  searchingZone[4];
	int workingPos[2];
	int prevWorkingPos[2];  /*
      for(i = 0; i < fl->nFeatures; i++)
	{
	  // fl->feature[i]->x += ob.workingPos[0];
	  // fl->feature[i]->y += ob.workingPos[1];
	  //  fl->feature[i]->errorx += ob.workingPos[0];
	  // fl->feature[i]->errory += ob.workingPos[1];
	  x   = (int) fl->feature[i]->errorx;
	  y   = (int) fl->feature[i]->errory;
	  // FIRST: Invert probability for all the points in its t-1 position (tracked or not)
	  _invertFeaturemap(x, y, featmap, tc->mindist, jImg2.width(), jImg2.height());
	  if (fl->feature[i]->val >= 0)  {
	    // SECOND: Update probability for all the tracked points in its current position
	    x   = (int) fl->feature[i]->x;
	    y   = (int) fl->feature[i]->y;
	    _fillFeaturemap(x, y, featmap, maskProb, tc->mindist, jImg2.width(), jImg2.height() );
	  }
	}

	 */

	int workingDims[2];
	int prevWorkingDims[2];

	//--- kalman prediction
	KalmanFilter *kalman; //! kalman filter (see filter jafar module)
	int estimatedPos[2]; //! position estimated by the kalman filter
	int estimatedSpeed[2]; //! speed estimated by the kalman filter
	int confidenceDist; //! confidence distance estimated by the kalman filter

	TrackedObject()
	{
		prevScope = 0;
		scope = 0;
		status = 1;
		flagRegion = 0;
		estimatedPos[0] = 0;
		estimatedPos[1] = 0;
		confidenceDist = 0;
		kalman = NULL;
		curFrame = 0;
		maxDisplacement = 30;  //50
		maxAccel = 0;  //50
		measureConfidence = 10;
		initSpeedConfidence = 10; //100
		securityMargin = 10;
		dt = 0.1;
	}

	TrackedObject(double datatime)
	{
		prevScope = 0;
		scope = 0;
		status = 1;
		flagRegion = 0;
		estimatedPos[0] = 0;
		estimatedPos[1] = 0;
		confidenceDist = 0;
		kalman = NULL;
		curFrame = 0;
		maxDisplacement = 30;  //50
		maxAccel = 20;  //50
		measureConfidence = 10;
		initSpeedConfidence = 10; //100
		securityMargin = 10;
		dt = datatime;

	}

	~TrackedObject()
	{
	}

	void writeKalmanToImg(cv::Mat &jRGBimg);
	void writeBordersToImg(cv::Mat &jRGBimg, KLT_TrackingContext tc);
	void writeContourToImg(cv::Mat &jRGBimg, unsigned char *featm);
	int detectLandscape (KLT_FeatureList feat);
	/**
	 @brief calculate the state of the object from the feature list (position, scope, boundingBox, etc)
	 @param fl the feature list
	 @param first specify if this is the initial selection
	 */
	int calcState();
	int calcConvex();
	/**
	 @brief detect and eliminate features which are not on the
	 object (with an heuristic on displacements and positions)
	 @param fl the feature list
	 */

};

typedef struct objM {
	tracker::TrackedObject * o;
	int numClus;
	int nP;
}objM, *objMptr;

typedef struct lstObj{
	objMptr ob;
	lstObj *prev;
	lstObj *next;
} lstObj, *lstObPtr;


typedef struct {
	lstObPtr tar;
	int totalClusters;
} lstclas;

typedef struct fusionObs{
	double **v;
	double state[4];
	double meanVel[2];
	int numClus;
	int nP;
}fusionObs, *fusObptr;

typedef struct lstFus {
	fusObptr ob;
	lstFus *prev;
	lstFus *next;
} lstFus, *lstFusptr;

int countNumberOfObjectPoints(mat _Xp,int ptsLastIm, int nc);
void verifyWorkingDims(int workingDims[2], int w, int h);
void verifyWorkingPos(int workingDims[2], int workingPos[2], int w, int h);
void orderNumberOfObjects(mat &_Xp, int nClus);
fusObptr fillStateObject (mat & obj, int nC);
int fusionTwoObjects(fusObptr ob1, fusObptr ob2);
double meanOfPoints(double **_Xpol, int nP, int coor);
double minimalPos(double **_Xpol, int nP, int coor);
double maximalPos(double **_Xpol, int nP, int coor);
lstFusptr findHeadFusionList(lstFusptr list);
lstObPtr findHeadOfList(lstObPtr list);
lstObPtr findTailOfList(lstObPtr list);
void shiftNumberOfObject(lstObPtr &list, int eraser);
void saveObjectsOfLinkedList(lstObPtr &lt, int nOb);
void addObjectsOfLinkedList(lstObPtr &lt, int nOb);
int extractPointsOfObject(mat &_Xp, mat & obj, int nc);
void eraseListObjects (lstObPtr & _eG);
int pointsToStructure(mat _Xp, fusObptr &aux, int nc);
void renameMatrixObjectIndex (mat &Xp, int name, int rename);
void deleteListObjects (lstFusptr & _eG);
void freeDoubleMatrix(double **m);
void initialiseStructObject(mat obj,  objMptr ob1, int numC, int nP, int ptsLastIm);
double meanOfOrientation(double **_Xpol, int nP);
int depthValidationOfGroups(mat & Xp, char *filename, int nClusters, int pActual);
} // end tracker

#endif /* MOT_HPP_ */
