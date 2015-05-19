/*
 * tracker.hpp
 *
 *  Created on: Feb 5, 2013
 *      Author: luzdora
 */

#ifndef TRACKER_HPP_
#define TRACKER_HPP_

#include "mot.hpp"
#include "klt.hpp"
//#include "regions.hpp"
extern "C" {
#include "klt.h"
}

#define MINI 1

  namespace tracker {
  using namespace boost_incl;

    int extractObjectData(mat _Xp, mat & obj,
			  int ptsLastIm, int nc, int ini, int shift);

    int jointObjectData(mat &_Xp, lstclas & ltc , mat &fus,
			int numClus, int totCl,  int ptsIm);

    int mergingObjects(mat & Xp, int nclus, int numClusBefore);

    void eraserlstclass (lstclas &ltc);

    int initTrackingObjects(lstclas &ltc, mat Xp,
			     vec & X0, mat & P0,
			     int nC, int nR, KLT_TrackingContext tc,
			     int ini, int numCls, int pts, double dt, int snkFlag);

    void initKalmanFilterObject(TrackedObject &ob, int tcx,
				int tcy, vec &xo,
				mat &po, int clus, int snkFlag);

    void kalmanFilterPrediction (lstclas &ltc, double dt,
				 int nC, int nR, KLT_TrackingContext tc,
				 cv::Mat &imR,
				 uchar *featuremap);
    void drawObjectDetection (lstclas &ltc, KLT_TrackingContext tc,
			      cv::Mat &imR,
			      uchar *featuremap);

    void kalmanPrediction (TrackedObject & ob, double dt,
			   int nC, int nR, KLT_TrackingContext tc );

    int objectTracking (lstclas &ltc, KLT_TrackingContext tc,
    		cv::Mat &jImg1, cv::Mat &jImg2,
			 vec & Z, mat H, vec & R,
			int nC, int nR, cv::Mat &imR,
			uchar * featuremap, uchar *maskProbability, double dt, int nIma);

    int modelDetection (TrackedObject & ob, KLT_TrackingContext tc,
    		cv::Mat &jImg1, cv::Mat &jImg2,
			uchar * featmap, uchar *maskProb, double dt);

    void updateKalmanFilter (TrackedObject & ob, vec & Z,
			     mat  H, vec & R );

    void deleteObjectZone (lstclas &ltc, uchar * featuremap, uchar * maskProb,
			   KLT_TrackingContext tc, int nC, int nR, int nIma);

    void writeKalmanObjectToImg(lstclas &ltc, cv::Mat &I,
				int first);

    void deleteObjectPoints (TrackedObject to, uchar * featuremap, uchar * maskProb,KLT_TrackingContext tc, int nC, int nR );

    void drawResults(TrackedObject & ob, cv::Mat &imR, uchar * featuremap, KLT_TrackingContext tc );

	int regionDetection(TrackedObject & ob, cv::Mat &jImg1,
						cv::Mat &jImg2, uchar *featmap,
			            uchar *maskProb, double dt, FILE *fp);

	void deleteObjectRegion (TrackedObject to, uchar * featuremap, int distMin, int nC, int nR );

    KLT_TrackingContext CreateTrackingContextAuxiliar(void);

    void  guardarDatosObjetos(int numob, int numIma, TrackedObject & ob);

  } // end tracker



#endif /* TRACKER_HPP_ */
