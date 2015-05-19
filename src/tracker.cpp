/*
 * tracker.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: luzdora
 */
#include "tracker.hpp"

namespace tracker {

using namespace klt;
using namespace tracker;
using namespace snake;
using namespace filter;
using namespace boost_incl;

int extractObjectData(mat _Xp, mat & obj, int ptsLastIm, int nc, int ini, int shift)
{
	int i, l, j;
	//  std::cout << "Size of xp: " << _Xp.size1() <<std::endl;
	for( i = _Xp.size1() - ptsLastIm, l=ini; i < (int)_Xp.size1(); i++)
	{
		if (_Xp(i,4) == nc)
		{
			for (j =0; j < (int)_Xp.size2(); j++)
				obj (l,j) = _Xp(i,j);

			obj (l,4)+=shift; // new index of object
			l++;
		}
	}
	//      std::cout << "Size of object extracted: " << l-ini <<std::endl;
	//  std::cout << "EXTRACT OBJ DATA Index: " << obj(ini,4) << " Inx " << obj(l-1,4) <<" Size obj " << l-ini << std::endl;
	return l;
}

// This function returns the number of objects fusionned (including the already existed clusters)
int jointObjectData(mat &_Xp, lstclas & ltc , mat &fus,
		int numClus, int totCl,  int ptsIm)
{
	int i, l, j, k, dim, nP=0, newCluster = 0;
	fus.clear(); // new fusion of points
	dim = _Xp.size2();

	// First extract points of the existant objects
	lstObPtr lt = ltc.tar;
	l=0;
	while (lt != NULL)
	{
		//  if ((*(lt->ob->o)).flagRegion != 1){
		for (i=0; i<lt->ob->nP; i++, l++){
			for (j=0; j< dim; j++) {
				fus(l,j) = (*(lt->ob->o)).fl(i,j);
				// std::cout << "Valor de fus " << fus(l,j) << std::endl;
			}
			fus(l,4) = lt->ob->numClus;
			//   std::cout << "Valor de fus " << fus(l,4) << std::endl;
		}
		// }
		lt = lt->next;
	}
	//   std::cout << "Size of fus object already detected " << l << std::endl;

	//  std::cout << "numClus " << numClus << " Total clus "<< totCl << std::endl;
	//second: Change the number of clusters
	j = totCl - numClus;
	for (k=1; k<=numClus; k++)
	{
		//	  std::cout << "Numero shift " << j << " para el cluster " << k << " punts plIm " << ptsIm << std::endl;
		nP = countNumberOfObjectPoints(_Xp, ptsIm, k);
		if ( nP > MINI)
		{
			fus.resize(l+nP,dim,1);
			l = extractObjectData(_Xp, fus, ptsIm, k, l, j);
			//    std::cout << "New size of objects already detected " << l << std::endl;
			newCluster++;
		}
		else
		{
			// One object less
			j--;
			std::cout<< "TRACKER JointData: There is not points in last image for object " << k << std::endl;
		}
	}

	if (newCluster==0)
		//	return (totCl - numClus); // not new objects
		return newCluster;
	else
		return (newCluster+(totCl - numClus));
}

void eraserlstclass ( lstclas &ltc)
{
	if(ltc.tar!=NULL){
	ltc.tar=findHeadOfList(ltc.tar);
	eraseListObjects (ltc.tar);
	}
	//   free (*ltc);
}

// This function returns the number of points of all the detected objects
int initTrackingObjects( lstclas & ltc, mat Xp, vec & X0, mat & P0,
		int nC, int nR, KLT_TrackingContext tc, int ini, int numCls,
		int pts, double dt, int snkFlag)
{
	lstObPtr ltob = ltc.tar;
	lstObPtr  aux_lt= NULL;
	int i, anP=0, nP=0, finC=numCls, elim=0;

	if (ini == 0) { // Initialization d'objects
		ltob = NULL;
		// Comenzamos el analisis para cada numero de cluster en Xp
		for (i=1; i<=numCls; i++)
		{
			anP = countNumberOfObjectPoints(Xp, pts, i);
			if (anP > MINI)
			{
				nP += anP; // add all the points that forms different objects
				//  std::cout << "Create objet : "<< i-elim << "with  points = " << anP << std::endl;
				aux_lt = (lstObPtr)calloc(1,sizeof(lstObPtr));
				aux_lt->ob = (objMptr)calloc(1,sizeof(objM));

				if (aux_lt->ob != NULL)
					aux_lt->ob->o = new TrackedObject(dt);
				else
				{
					std::cout << "Problems, No memory reserved for objMptr" << std::endl;
					return 0;
				}

				initialiseStructObject(Xp, aux_lt->ob, i, anP, pts);
				aux_lt->ob->numClus = i-elim;
				if(ltob==NULL)
				{
					aux_lt->prev = NULL;
					aux_lt->next = NULL;
					ltob = aux_lt;
				}
				else
				{
					ltob->next = aux_lt;   //creo que esta linea esta de mas
					aux_lt->prev = ltob;
					ltob = aux_lt;
					ltob->next = NULL;
				}
				initKalmanFilterObject(*(ltob->ob->o), tc->borderx, tc->bordery, X0, P0,i-elim, snkFlag);
				//	  verifyWorkingDims((*(ltob->ob->o)).workingDims, nC, nR);
			}
			else
			{
				if (numCls == 1)
					return 0;
				elim++;
				finC--; //eliminate one cluster because it doesnt containt points
				std::cout << " INIT TRACKING: There is not pts in last image for cluster " << i << "Total Clusters " << finC << std::endl;
				if (finC == 0)
					return 0;
			}
		}//end for
		ltob = findHeadOfList(ltob);
		ltc.tar = ltob;
		ltc.totalClusters = finC;
	}
	else  // One or more objects were already tracked
	{
		orderNumberOfObjects(Xp, numCls);
		// take the head of the list
		ltob = findHeadOfList(ltob);
		aux_lt = ltob;  // save a copy
		if (numCls < ini)
			saveObjectsOfLinkedList(aux_lt, numCls);
		else if (numCls > ini)
			addObjectsOfLinkedList(aux_lt, numCls-ini);
		// At this point, the list has the same number of nodes as objects
		aux_lt = ltob; // head of list
		i=0;
		while (aux_lt != NULL)
		{
			i++;
			anP = countNumberOfObjectPoints(Xp, pts, i);
			nP += anP; // add all the points that forms different objects
			initialiseStructObject(Xp, aux_lt->ob, i, anP, pts);
			aux_lt->ob->numClus = i;
			initKalmanFilterObject(*(aux_lt->ob->o), tc->borderx, tc->bordery, X0, P0,i, snkFlag);
			//   verifyWorkingDims((*(aux_lt->ob->o)).workingDims, nC, nR);
			aux_lt = aux_lt->next;
		}
		std::cout << "The number of objects save " << i << std::endl;
		ltob = findHeadOfList(ltob);
		ltc.tar = ltob;
		ltc.totalClusters = i;
	}
	// The total number of points taked to form the different objects
	//     std::cout << "Bye Bye nP = " << nP << std::endl;
	return nP;
}

void initKalmanFilterObject(TrackedObject &ob, int tcx, int tcy, vec &xo,
		mat &po, int clus, int snkFlag)
{
	ob.status = 1;
	ob.flagRegion = snkFlag;
	int points = ob.calcState(); // initialization

	/****** SNAKE CONTOUR : HAY QUE ACTIVAR BANDERA    ******/

	if (ob.flagRegion && points>4 ){
		//if (ob.conv.size1()>4 ){
		if(ob.calcConvex()>3){
			snake::findContour(ob.conv, ob.contour);
			ob.flagRegion = 1; // Hay snake
		} else  {
			std::cout << " INITIALIZATION : Initial Contour will be the bounding-Box " << std::endl;
			ob.contour.resize(4,2);
			ob.contour(0,0)= ob.boundingBox[x1p];
			ob.contour(0,1)= ob.boundingBox[y1p];
			ob.contour(1,0)= ob.boundingBox[x2p];
			ob.contour(1,1)= ob.boundingBox[y2p];
			ob.contour(2,0)= ob.boundingBox[x2p];
			ob.contour(2,1)= ob.boundingBox[y1p];
			ob.contour(3,0)= ob.boundingBox[x1p];
			ob.contour(3,1)= ob.boundingBox[y2p];
			//	ob.flagRegion = 0;  // No contour snake
			//		ob.conv.resize(0,0,0);
			//		ob.contour.resize(0,0,0);
		}
	} else
		ob.flagRegion = 0;
	// Buscar la zona rectangular del objeto sea snake o no
	ob.workingDims[0] = ob.boundingBox[x2p] - ob.boundingBox[x1p] + 1 + ob.securityMargin + ob.maxDisplacement + 2*tcx;
	ob.workingDims[1] = ob.boundingBox[y2p] - ob.boundingBox[y1p] + 1 + ob.securityMargin + ob.maxDisplacement + 2*tcy;
	ob.workingPos[0]  = ob.pos[0]-ob.workingDims[0]/2;
	ob.workingPos[1]  = ob.pos[1]-ob.workingDims[1]/2;
	// Vector de estado
	//std::cout<<"  PrevWorkingPosX  "<<ob.prevWorkingPos[0]<<"  PrevWorkingPosY  "<<ob.prevWorkingPos[1]<<"  WorkingPosX  "<<ob.workingPos[0]<<"  WorkingPosY  "<<ob.workingPos[1]<<std::endl;
	xo(0) = ob.pos[0]; xo(1) = ob.pos[1]; xo(2) =ob.vel[0] ; xo(3) = ob.vel[1];
	//double p = ob.measureConfidence*ob.measureConfidence;
	//    double p = ob.measureConfidence;
	//   std::cout << "Desde Init Kalman Object Vel X " << ob.vel[0] << " VY " << ob.vel[1] << " Of object " << clus << std::endl;
	//double vx = ob.vel[0]*ob.vel[0];
	//double vy = ob.vel[1]*ob.vel[1];

	//double vx = 5*5;
	//double vy = 5*5;
	//  double v = ob.initSpeedConfidence*ob.initSpeedConfidence;
	po(0,0) = 1; po(0,1)=0; po(0,2)=0; po(0,3)=0;
	po(1,0) = 0; po(1,1)=1; po(1,2)=0; po(1,3)=0;
	po(2,0) = 0; po(2,1)=0; po(2,2)=1; po(2,3)=0;
	po(3,0) = 0; po(3,1)=0; po(3,2)=0; po(3,3)=1;

	ob.kalman = new KalmanFilter(4);
	ob.kalman->init(xo, po);
}

void kalmanFilterPrediction (lstclas &ltc, double dt, int nC, int nR, KLT_TrackingContext tc, cv::Mat &imR, uchar *featuremap)
{
	lstObPtr lt = ltc.tar;
	while (lt != NULL)
	{
		kalmanPrediction(*(lt->ob->o), dt, nC, nR, tc );
		 drawResults(*(lt->ob->o), imR, featuremap, tc);
		lt = lt->next;
	}
}

void kalmanPrediction (TrackedObject & ob, double dt, int nC, int nR, KLT_TrackingContext tc)
{
	// Dims y pos  ahora anteriores
	ob.prevWorkingDims[0] = ob.workingDims[0];
	ob.prevWorkingDims[1] = ob.workingDims[1];
	ob.prevWorkingPos[0]  = ob.workingPos[0];
	ob.prevWorkingPos[1]  = ob.workingPos[1];

	mat F(4,4);
	F(0,0)=1; F(0,1)=0; F(0,2)=dt; F(0,3)=0;
	F(1,0)=0; F(1,1)=1; F(1,2)=0;  F(1,3)=dt;
	F(2,0)=0; F(2,1)=0; F(2,2)=1;  F(2,3)=0;
	F(3,0)=0; F(3,1)=0; F(3,2)=0;  F(3,3)=1;


	// boost::numeric::ublas::sym_mat<double> Q(4,4);
	sym_mat Q(4,4);
	//double a = ob.maxAccel*ob.maxAccel;

	double acel[] = {(ob.vel[0]- ob.preVel[0])/dt,(ob.vel[1]- ob.preVel[1])/dt};


	double a = acel[0]*acel[0]+acel[1]*acel[1];
	// error of the model = acceleration a != 0 => error on speed = a.dt => error on pos = Sa.dt=a.dt^2/2
	Q(0,0)=a*pow(dt,4)/4;	Q(0,1)=0;              	Q(0,2)=a*pow(dt,3)/2; 	Q(0,3)=0; // delta t = 1 -> pixels/interframe
	Q(1,0)=0;		        Q(1,1)=a*pow(dt,4)/4; 	Q(1,2)=0;             	Q(1,3)=a*pow(dt,3)/2;
	Q(2,0)=a*pow(dt,3)/2;	Q(2,1)=0;              	Q(2,2)=a*pow(dt,2);   	Q(2,3)=0;
	Q(3,0)=0;    	        Q(3,1)=a*pow(dt,3)/2; 	Q(3,2)=0;             	Q(3,3)=a*pow(dt,2);

	LinearPredictModel predictModel(F, Q);
	ob.kalman->predict(predictModel);

	// ******************** use kalman prediction
	vec &X = ob.kalman->getX();


	ob.estimatedPos[0]   = (int)X(0); ob.estimatedPos[1]   = (int)X(1);
	ob.estimatedSpeed[0] = (int)X(2); ob.estimatedSpeed[1] = (int)X(3);



	//jblas::sym_mat &P = ob.kalman->getP();
	// boost::numeric::ublas::sym_mat<double> &P = ob.kalman->getP();
	sym_mat &P = ob.kalman->getP();

	if (P(0,0) > P(1,1)) ob.confidenceDist = (int)sqrt(P(0,0)); else ob.confidenceDist = (int)sqrt(P(1,1));

	//     std::cout << "In prediccion confidence dist : "<< ob.confidenceDist  << std::endl;
	//    std::cout << " EstimatedSpeed : "<<ob.estimatedSpeed[0]<<", "<<ob.estimatedSpeed[1] << std::endl;

	// calc working pos de acuerdo con la prediction

	ob.workingPos[0] = ob.estimatedPos[0]-ob.workingDims[0]/2;
	ob.workingPos[1] = ob.estimatedPos[1]-ob.workingDims[1]/2;
	verifyWorkingPos(ob.workingDims, ob.workingPos, nC, nR);

	//    std::cout << "WorkingPos : " << ob.workingPos[0] <<", "<< ob.workingPos[1] <<", EstimatedPos : "<<ob.estimatedPos[0]<<", "<<ob.estimatedPos[1] << std::endl;
	//std::cout << "PrevWorkingPos : " << ob.prevWorkingPos[0] <<", "<< ob.prevWorkingPos[1] << std::endl;
	//  std::cout << "WorkingDims : " << ob.workingDims[0] <<", "<< ob.workingDims[1] <<", PrevWorkDims : "<<ob.prevWorkingDims[0]<<", "<<ob.prevWorkingDims[1] << std::endl;

}

void drawObjectDetection (lstclas &ltc, KLT_TrackingContext tc, cv::Mat &imR, uchar *featuremap)
{
	lstObPtr lt = ltc.tar;
	while (lt != NULL)
	{
		drawResults(*(lt->ob->o), imR, featuremap, tc);
		lt = lt->next;
	}
}

void drawResults(TrackedObject & ob, cv::Mat &imR, uchar * featuremap, KLT_TrackingContext tc)
{
	if (ob.status){
		ob.writeKalmanToImg(imR); // first time
		ob.status = 0;
	}
	else  {
		//if ( (ob.conv.size1()>4)  && (ob.conv.size1()<10) && (ob.flagRegion)){
		if ( (ob.conv.size1()>4) && (ob.flagRegion)){
			std::cout << "Deja a Matlab que dibuje el contorno...jejejeje" << std::endl;
			/*	 ob.calcConvex();
   	 // if (ob.conv.size1()>4 && (ob.conv.size1()<10)){
   	  if (ob.conv.size1()>4){
   	    findContour(ob.conv, ob.contour);
   	  // At least four points are needed for finding a contour
   	  ob.writeContourToImg(imR, featuremap);
   	   AQUI FALTAN LLAVES...OJO
   	  }
   	  else
   	    ob.writeBordersToImg(imR, tc);
   	} */
		}
		else
			ob.writeBordersToImg(imR, tc);

	}
}

void guardarDatosObjetos(int numob, int numIma, TrackedObject & ob)
{
	FILE *fp;
	int i;

	fp = fopen("/tmp/objects.txt", "a");

	for( i = 0; i < (int)ob.fl.size1(); i++)
		fprintf(fp,"%d\t%d\t%5.1f\t%5.1f\t%5.1f\t%5.1f\n", numIma, numob, (float) ob.fl(i,0),(float)ob.fl(i,1), (float)ob.fl(i,2), (float)ob.fl(i,3));

	fclose(fp);

}


int objectTracking (lstclas & ltc, KLT_TrackingContext tc,
		cv::Mat &jImg1, cv::Mat &jImg2,
		vec & Z, mat H, vec & R,
		int nC, int nR, cv::Mat &imR,
		uchar *featuremap, uchar *maskProbability, double dt, int nIma)
{
	lstObPtr lt = ltc.tar;
	int ntrakP=0, shif=0;

	FILE *fp=NULL;
	double t = 100000;
	char fname1[100] = "/tmp/banderac.txt";
	char fname2[100] = "/tmp/snake.txt";
	int numob = 0;

	while (lt != NULL && fp == NULL)
	{
		numob++;
		guardarDatosObjetos(numob, nIma, *(lt->ob->o));

		if (lt->ob->o->flagRegion) {
			while ((fp = fopen(fname1 , "rb")) != NULL){
				fclose(fp);
				usleep(t);
			}
			fp = fopen(fname2 , "r");
			if (fp == NULL){
				printf("Adiosin...No hay nada huesin...");
				exit(0);
			}
		}
		lt = lt->next;
	}
	lt = ltc.tar;

	while (lt != NULL)
	{
		if (shif>0)
			lt->ob->numClus -=shif;

		if (lt->ob->o->flagRegion)
			lt->ob->nP = regionDetection(*(lt->ob->o), jImg1, jImg2, featuremap, maskProbability, dt, fp);
		else
			lt->ob->nP = modelDetection(*(lt->ob->o), tc, jImg1, jImg2, featuremap, maskProbability, dt);
		//  lt->ob->nP = regionDetection(*(lt->ob->o), jImg1, jImg2, featuremap, maskProbability, dt, "/tmp/mask.bin");

		if (lt->ob->nP <= MINI) // Groups need to have more than a certain number of points
		{
			lstObPtr aux;
			if (ltc.totalClusters >= 2) {
				// Verify if there is the head
				if (lt->prev == NULL)
				{
					aux = lt->next;
					lt->next = NULL;
					eraseListObjects (lt);
					aux->prev = NULL; //new head
					ltc.tar = aux;
					lt = aux;
				}
				else if (lt->next == NULL) // if there is the tail
				{
					aux = lt->prev;
					lt->prev = NULL;
					eraseListObjects (lt);
					aux->next = NULL;
					lt = aux->next;
				}
				else{
					/*************  AQUI ME QUEDE FALTA ELIMINAR EL OBJETO QUE ESTA EN MEDIO DE LA LISTA Y ASEGURAR QUE LA LISTA QUEDE BIEN ENLAZADA ***/
					//   std::cout << "More than 2 objects " << std::endl;
					//Unlink the element
					aux = lt->prev;
					lt->prev = NULL;
					aux->next = lt->next;
					lt->next = NULL;
					aux->next->prev = aux;
					eraseListObjects (lt);
					lt = aux->next;
				}
				//  std::cout << "Adios al primero " << std::endl;
				ltc.totalClusters -= 1;
				shif++;
				continue;
			}
			else
			{
				std::cout << " No more objects for tracking in the scene " << std::endl;
				eraseListObjects (lt);
				ltc.tar = NULL;
				ltc.totalClusters = 0;
				return 0;
			}
		}
		else {  // Group with a minimal number of points
			ntrakP += lt->ob->nP;
			updateKalmanFilter(*(lt->ob->o), Z, H, R);
			kalmanPrediction(*(lt->ob->o), dt, nC, nR, tc);
 			drawResults(*(lt->ob->o), imR, featuremap, tc);

			lt = lt->next;
		}
	}
	if (fp!=NULL)
		fclose(fp);
	return ntrakP;
}


int regionDetection(TrackedObject & ob, cv::Mat &jImg1,
		cv::Mat &jImg2, uchar *featmap,
		uchar *maskProb, double dt, FILE *fp)
{

	fscanf(fp, "%lf %lf %lf %lf %d %d %d %d\n", &ob.pos[0], &ob.pos[1],
			&ob.vel[0], &ob.vel[1], &ob.boundingBox[0], &ob.boundingBox[1],
			&ob.boundingBox[2], &ob.boundingBox[3]);

	// Esto no creo que sea necesrio xq el Matlab se encarga de to_do, solo para que no se queje
	ob.workingDims[0] = ob.boundingBox[x2p] - ob.boundingBox[x1p] + 1 + ob.securityMargin + ob.maxDisplacement ;
	ob.workingDims[1] = ob.boundingBox[y2p] - ob.boundingBox[y1p] + 1 + ob.securityMargin + ob.maxDisplacement ;
	ob.workingPos[0]  = ob.pos[0]-ob.workingDims[0]/2;
	ob.workingPos[1]  = ob.pos[1]-ob.workingDims[1]/2;

	ob.searchingZone[0] = (int) ob.boundingBox[0]-10;
	ob.searchingZone[1] = (int) ob.boundingBox[1]-10;
	ob.searchingZone[2] = (int) ob.boundingBox[2]+10;
	ob.searchingZone[3] = (int) ob.boundingBox[3]+10;

	return (int) ob.fl.size1();

}

//// MODIFICAME Y RAPIDO!!!!!
int modelDetection (TrackedObject & ob, KLT_TrackingContext tc,
		cv::Mat &jImg1, cv::Mat &jImg2,
		uchar *featmap, uchar *maskProb, double dt)
{
	KLT_FeatureList fl;
	int i, j,nFeat, nLostFeatures=0;
	int *index = NULL;
	fl=KLTCreateFeatureList((int)ob.fl.size1());

	// Verification du WorkingPos
	if ( ob.prevWorkingPos[0] < ob.workingPos[0]-10 ){
		std::cout << "PrevWorkingPosX : " << ob.prevWorkingPos[0] << std::endl;
		ob.prevWorkingPos[0] = ob.workingPos[0];
		std::cout << "NOW!! PrevWorkingPosX : " << ob.prevWorkingPos[0] << std::endl;
	}
	if ( ob.prevWorkingPos[1] < ob.workingPos[1]-10 ){
		std::cout << "PrevWorkingPosY : " << ob.prevWorkingPos[1] << std::endl;
		ob.prevWorkingPos[1] = ob.workingPos[1];
		std::cout << "NOW!! PrevWorkingPosY : " << ob.prevWorkingPos[1] << std::endl;
	}

    //-------------- Verify image dimensions
	ob.prevWorkingPos[1] = ob.prevWorkingPos[1] < 0 ? 0 : ob.prevWorkingPos[1];
	ob.prevWorkingPos[0] = ob.prevWorkingPos[0] < 0 ? 0 : ob.prevWorkingPos[0];

	ob.workingPos[0] = ob.workingPos[0] < 0 ? 0 : ob.workingPos[0];
	ob.workingPos[1] = ob.workingPos[1] < 0 ? 0 : ob.workingPos[1];

	int maxy = max(ob.workingPos[1],ob.prevWorkingPos[1]);
	int maxx = max(ob.workingPos[0],ob.prevWorkingPos[0]);

	if(maxy + ob.workingDims[1] > jImg1.rows)
		ob.workingDims[1] = jImg1.rows - maxy;

	if(maxx + ob.workingDims[0] > jImg1.cols)
		ob.workingDims[0] = jImg1.cols - maxx;


	// ******************** adjust coordinates
	for( i = 0; i < (int)ob.fl.size1(); i++)
	{
		fl->feature[i]->x = ob.fl(i,0);
		fl->feature[i]->y = ob.fl(i,1);
		fl->feature[i]->val = ob.fl(i,4);
		fl->feature[i]->auxval = ob.fl(i,5);
		// change coordinate references
		fl->feature[i]->x -= ob.prevWorkingPos[0];
		fl->feature[i]->y -= ob.prevWorkingPos[1];
	}

	// ******************** get subimage where to search, with conversion between jafar image and aligned data
	KLT_PixelType *img1_ = NULL, *img2_ = NULL;
	cv::Mat img1, img2;

	// if (tc->pyramid_last == NULL)

	// img1 = compatImage(jImg1, treatment1, ob.prevWorkingPos[0], ob.prevWorkingPos[1], ob.workingDims[0], ob.workingDims[1]); // (verification required ...)
	// img2 = compatImage(jImg2, treatment2, ob.workingPos[0], ob.workingPos[1], ob.workingDims[0], ob.workingDims[1]);


	img1 = jImg1(cv::Range(ob.prevWorkingPos[1], ob.prevWorkingPos[1]+ob.workingDims[1]),
    		     cv::Range(ob.prevWorkingPos[0], ob.prevWorkingPos[0]+ob.workingDims[0]));

    img2 = jImg2(cv::Range(ob.workingPos[1], ob.workingPos[1]+ob.workingDims[1]),
    		     cv::Range(ob.workingPos[0], ob.workingPos[0]+ob.workingDims[0]));



	img1_ = img1.data;
	img2_ = img2.data;

	//std::cout<<"  PrevWorkingPosX  "<<ob.prevWorkingPos[0]<<"  PrevWorkingPosY  "<<ob.prevWorkingPos[1]<<"  WorkingPosX  "<<ob.workingPos[0]<<"  WorkingPosY  "<<ob.workingPos[1]<<std::endl;

	KLTTrackFeatures(tc, img1_, img2_, ob.workingDims[0], ob.workingDims[1], fl);

	// ******************** calculate object.boundingBox and object.pos
	// change coordinate references for features and adjust bounding box

	for(i = 0; i < fl->nFeatures; i++)
	{
		fl->feature[i]->x += ob.workingPos[0];
		fl->feature[i]->y += ob.workingPos[1];
		fl->feature[i]->errorx += ob.workingPos[0];
		fl->feature[i]->errory += ob.workingPos[1];
	}
	int lost = ob.detectLandscape(fl);
	nLostFeatures = fl->nFeatures - KLTCountRemainingFeatures(fl);

	if (nLostFeatures > 0 ){
		int flag, k=0, dim = ob.fl.size2();
		mat aux (fl->nFeatures, dim);
		aux = ob.fl;
		ob.fl.resize(fl->nFeatures-nLostFeatures, dim, 0);
		index = (int *) KLTFindIndexOfLostFeatures(fl, nLostFeatures);
		for(i = 0; i < fl->nFeatures; i++)
		{
			flag = 1;
			for (j=0; j< nLostFeatures; j++)
				if (i == *(index+j)){
					flag = 0;
					break;
				}
			if (flag){
				for (j=0; j < dim; j++)
					ob.fl(k,j) = aux(i,j);
				k++;
			}
		}
		free(index);
	}
	// If there are any point in the object set
	if (lost >= 0){
		ob.curFrame++;
		int k=0;
		// ******************** calculate object.boundingBox and object.pos
		// change coordinate references for features and adjust bounding box
		for(i = 0, k=0; i < fl->nFeatures; i++)
		{
			// antes calculo la velocidad
			if (fl->feature[i]->val == 0){
				ob.fl(k,2) =  (fl->feature[i]->x - ob.fl(k,0))/dt;
				ob.fl(k,3) =  (fl->feature[i]->y - ob.fl(k,1))/dt;
				// std::cout << " ob.fl(2) " << ob.fl(k,2) << " DE " << fl->feature[i]->x << " Y " <<  ob.fl(k,0) << std::endl;
				ob.fl(k,0) = fl->feature[i]->x;
				ob.fl(k,1) = fl->feature[i]->y;
				ob.fl(k,4) = fl->feature[i]->val;
				//      std::cout << " points[" << i << "].x = " <<  ob.fl(k,0) << " ; " << std::endl;
				//  std::cout << " points[" << i << "].y = "  << ob.fl(k,1) << " ; " << std::endl;
				k++;
			}
			else if (fl->feature[i]->val > 0)
			{
				ob.fl(k,0) = fl->feature[i]->x;
				ob.fl(k,1) = fl->feature[i]->y;
				ob.fl(k,2) =  0.0;
				ob.fl(k,3) =  0.0;
				ob.fl(k,4) = fl->feature[i]->val;
				ob.fl(k,5) = fl->feature[i]->auxval;
				k++;
			}
		}

		nFeat = ob.calcState();

		//RE-initializacion de las dimensiones de busqueda
		// the dimensions of the working area
		/* SI REGIONES  if (ob.flagRegion && ob.conv.size1()>4)
	  {
	    ob.workingDims[0] = ob.boundingBox[x2p] - ob.boundingBox[x1p] + 1 + 2*ob.securityMargin;
	    ob.workingDims[1] = ob.boundingBox[y2p] - ob.boundingBox[y1p] + 1 + 2*ob.securityMargin;
	    ob.workingPos[0]  = ob.pos[0]-ob.workingDims[0]/2;
	    ob.workingPos[1]  = ob.pos[1]-ob.workingDims[1]/2;
	    std::cout << " Set Contour desde model detection, numFeatures= " << nFeat << std::endl;
	    ob.calcConvex();

	    if (ob.conv.size1() > 4 ){
	    findContour(ob.conv, ob.contour);
	    setContourAndRegionInMask( ob.contour,ob.workingDims[0], ob.workingDims[1],
	    ob.workingPos[0], ob.workingPos[1], featmap, jImg2.width(), jImg2.height() );
	    }
	    else
	    {
	    std::cout << " Reset the size of conv y contour " << std::endl;
	    ob.conv.resize(0,0,0);
	    ob.contour.resize(0,0,0);
	    }

	  }
	  else{
	  std::cout << " Delete contour zone & reset the size of conv y contour " << std::endl;
	  ob.conv.resize(0,0,0);
	  ob.contour.resize(0,0,0);
		 */
		ob.workingDims[0] = ob.boundingBox[x2p] - ob.boundingBox[x1p] + 1 + ob.securityMargin + ob.maxDisplacement + 2*tc->borderx;
		ob.workingDims[1] = ob.boundingBox[y2p] - ob.boundingBox[y1p] + 1 + ob.securityMargin + ob.maxDisplacement + 2*tc->bordery;
		ob.workingPos[0]  = ob.pos[0]-ob.workingDims[0]/2;
		ob.workingPos[1]  = ob.pos[1]-ob.workingDims[1]/2;
		/*	}  SI REGIONES */

		ob.searchingZone[0] = (int) ob.boundingBox[0];
		ob.searchingZone[1] = (int) ob.boundingBox[1];
		ob.searchingZone[2] = (int) ob.boundingBox[2];
		ob.searchingZone[3] = (int) ob.boundingBox[3];

		/* ESTO DEBERIA IR EN LA DETECTION DE REGION */

		/*	if (nFeat>4 ){
	  //    if(){
	  std::cout << " Set Contour desde model detection, numFeatures= " << nFeat << std::endl;
	  ob.calcConvex();
	  if (ob.conv.size1() > 4 ){
	    findContour(ob.conv, ob.contour);
	    setContourAndRegionInMask( ob.contour,ob.workingDims[0], ob.workingDims[1],
				       ob.workingPos[0], ob.workingPos[1], featmap, jImg2.width(), jImg2.height() );
	  }
	  else
	    {
	      std::cout << " Reset the size of conv y contour " << std::endl;
	      ob.conv.resize(0,0,0);
	      ob.contour.resize(0,0,0);
	    }
	}
	else  // Delete contour zone {
	    std::cout << " Delete contour zone & reset the size of conv y contour " << std::endl;
	    ob.conv.resize(0,0,0);
	    ob.contour.resize(0,0,0);
	  }
		 */

	}
	else
	{
		std::cout << "All features in this object were lost " << std::endl;
		nFeat = 0;
	}

	KLTFreeFeatureList(fl);
	// destroyCompatImage(img1, treatment1);
	//   if (nLostFeatures == 0)
	// destroyCompatImage(img2, treatment2);
	return nFeat;
}

void updateKalmanFilter (TrackedObject & ob, vec & Z, mat H, vec & R )
{
	// ******************** update kalman
	//   jblas::vec R(2);
	//    R(0)=R(1)=ob.measureConfidence*ob.measureConfidence; // pixel precision for observation
	R(0)=R(1)=1; // pixel precision for observation
	LinearObserveModel observeModel(H);
	observeModel.setUncorrelatedR(R);
	//  jblas::vec Z(2);
	Z(0) = ob.pos[0];
	Z(1) = ob.pos[1];
	ob.kalman->update(observeModel, Z);
}

void deleteObjectZone (lstclas &ltc, uchar * featuremap, uchar *maskProb,
		KLT_TrackingContext tc, int nC, int nR, int nIma)
{
	lstObPtr lt = ltc.tar;
	FILE *fp, *fp2;
	char fname[40] = "/tmp/datosOb.txt";
	char fname1[40] = "/tmp/banderac.txt";
	char fname2[40] = "/tmp/dinap.txt";
	int ind=0;
	// Abrir archivo datosOb
	fp = fopen(fname , "w");
	if (fp == NULL){
		printf(" SNAKE: Can't open file '%s' for writing", fname);
		exit(0);
	}
	fprintf(fp, "%d\n", nIma);

	// Abrir archivo dinap
	fp2 = fopen(fname2 , "w");
	if (fp2 == NULL){
		printf(" SNAKE: Can't open file '%s' for writing", fname2);
		exit(0);
	}

	while (lt != NULL)
	{
		if ((*(lt->ob->o)).flagRegion == 1)
		{
			ind++;
			// Imagen completa
			if ((*(lt->ob->o)).status) //{ Primera vez, aun no hay contorno Matlab
				setContourAndRegionInMask((*(lt->ob->o)).contour, (*(lt->ob->o)).workingDims[0], (*(lt->ob->o)).workingDims[1], (*(lt->ob->o)).workingPos[0], (*(lt->ob->o)).workingPos[1], featuremap, nC, nR, ind);
			//	snake::setContourAndRegionInMask((*(lt->ob->o)).contour,(*(lt->ob->o)).workingDims[0], (*(lt->ob->o)).workingDims[1], (*(lt->ob->o)).workingPos[0], (*(lt->ob->o)).workingPos[1], featuremap, nC, nR, ind);
			//			std::cout<< " Primera vez, aun no hay contorno Matlab" << std::endl;
			else
				deleteObjectRegion(*(lt->ob->o), featuremap, tc->mindist, nC, nR );

			fprintf(fp, "%d %d %d %d\n", (*(lt->ob->o)).workingPos[0], (*(lt->ob->o)).workingPos[1],
					(*(lt->ob->o)).workingDims[0], (*(lt->ob->o)).workingDims[1]);
			fprintf(fp, "%lf %lf \n", (*(lt->ob->o)).estimatedPos[0]-(*(lt->ob->o)).pos[0], (*(lt->ob->o)).estimatedPos[1]-(*(lt->ob->o)).pos[0]);

			//Solo una region donde esta el objeto
			//  setContourAndRegionInMask((*(lt->ob->o)).contour,(*(lt->ob->o)).workingDims[0], (*(lt->ob->o)).workingDims[1],
			// 				0, 0, featuremap, (*(lt->ob->o)).workingDims[0], (*(lt->ob->o)).workingDims[1]);
			//      printf ("\n In contour Region ");
		}
		/*  DESCOMENTAR ESTO SI SE REQUIEREN LOS OBJETOS CON PUNTOS SIN CONTORNOS
	    else {

	      fprintf(fp2, "%d\n", (int)(*(lt->ob->o)).fl.size1() );
	      for (int i=0; i < (int)(*(lt->ob->o)).fl.size1() ; i++){
	    	x1 = 	( ( (*(lt->ob->o)).fl(i,0) ) + 0.5);
	    	y1 =   ( ( (*(lt->ob->o)).fl(i,1) ) + 0.5);
		fprintf(fp2, " %d %d ",  x1, y1 ); //Puntos del objeto no region
	      }
	      fprintf(fp2, "\n" ); //Puntos del objeto no region
	      deleteObjectPoints(*(lt->ob->o), featuremap, maskProb, tc, nC, nR);
	      }  */
		//	   std::cout<< " Delete zone of cluster: " << lt->ob->numClus << std::endl;
		lt = lt->next;
	}

	fclose(fp);
	fclose(fp2);

	// Genera archivo bandera indicador a matlab
	fp = fopen(fname1 , "w");
	if (fp == NULL){
		printf(" SNAKE: Can't open file '%s' for flag", fname);
		exit(0);
	}
	fclose(fp);
}

void writeKalmanObjectToImg(lstclas &ltc, cv::Mat &I, int first)
{
	lstObPtr lt = ltc.tar;
	while (lt != NULL)
	{
		(*(lt->ob->o)).writeKalmanToImg(I);
		//	   std::cout<< " Writing image zone of cluster: " << lt->ob->numClus << std::endl;
		lt = lt->next;
	}
}

void deleteObjectPoints (TrackedObject to, uchar * featuremap, uchar *maskProb,  KLT_TrackingContext tc, int nC, int nR )
{
	int j;
	int mind;
	double x1, y1;

	mind = tc->mindist;
	//  mind = 20;

	x1 = to.fl(0,0);
	y1 = to.fl(0,1);

	//otra version para evaluar los limites de la zona del objeto
	if (x1-mind<0)
		x1=0;
	else if (x1+mind >= nC)
		x1 = nC - mind;

	if (y1-mind<0)
		y1=0;
	else if (y1+mind >= nR)
		y1 = nR - mind;
	//    std::cout << " x1 " << x1 << " y1 " << y1 << " x2  " << x1+scopX << " y2 " << y1+scopY << std::endl;

	for(j=0; j < (int)to.fl.size1(); j++){
		x1 = to.fl(j,0);
		y1 = to.fl(j,1);
		_fillFeaturemap(x1, y1, featuremap, maskProb, mind, nC, nR);
	}
	// std::cout << " x1 = " << fl->feature[n]->x << " y1 " <<  fl->feature[n]->y << "num points n = " << n-1 << std::endl;
}

/*
    void drawResultsMatlab(TrackedObject & ob, image::Image &imR, uchar * featuremap, KLT_TrackingContext tc)
       {
         if (ob.status){
           ob.writeKalmanToImg(imR); // not first time
           ob.status = 0;
         }
         else  {
   	if ( (ob.conv.size1()>4)  && (ob.conv.size1()<10) && (ob.flagRegion)){
   	  ob.calcConvex();
   	  if (ob.conv.size1()>4 && (ob.conv.size1()<10)){
   	    findContour(ob.conv, ob.contour);
   	  // At least four points are needed for finding a contour
   	  ob.writeContourToImg(imR, featuremap);
   	  }
   	  else
   	    ob.writeBordersToImg(imR, tc);
   	}
   	else
   	  ob.writeBordersToImg(imR, tc);
         }
       }
 */
// Si no hay intercambio con Matlab

void deleteObjectRegion (TrackedObject to, uchar * featuremap, int distMin, int nC, int nR )
{
	int i, j;
	int mind;
	double x1, y1;

	double scopX = to.boundingBox[2]-to.boundingBox[0] + 1;
	double scopY = to.boundingBox[3]-to.boundingBox[1] + 1;

	mind = distMin;
	mind--;
	//  scop = to.scope;
	x1 = to.pos[0]-(scopX/2);
	y1 = to.pos[1]-(scopY/2);

	//otra version para evaluar los limites de la zona del objeto
	if (x1<0)
		x1=0;
	else if (x1+scopX >= nC)
		x1 = nC - scopX;

	if (y1<0)
		y1=0;
	else if (y1+scopY >= nR)
		y1 = nR - scopY;

	//     std::cout << " x1 " << x1 << " y1 " << y1 << " x2  " << x1+scopX << " nC " << nC << "nR" << nR << std::endl;
	/*
      for(j=0; j <=(int)(scopX/2*mind); j++)
	for (i=0; i<=(int)(scopY/2*mind) ; i++)
	  _fillFeaturemap(x1 + (2*mind*j), y1 + (2*mind*i), featuremap, maskProb, mind, nC, nR);
	 */
	//int corner = y1*nC + x1;
	for(i=x1; i <=(int)x1+scopX; i++)
		for (j=y1; j<=(int) y1+scopY; j++)
			featuremap[j*nC + i ] = 255;

	//	    std::cout << " x1 " << fl->feature[n]->x << " y1 " <<  fl->feature[n]->y << "num points n = " << n-1 << std::endl;

}


/* Funcion auxiliar para que funcione el trackingContext  */
KLT_TrackingContext CreateTrackingContextAuxiliar(void)
{
	return KLTCreateTrackingContext();
}

// This function returns the number of objects finally identified including detected objects,
// otherwise it returns 0 added objects.

int mergingObjects(mat & Xp, int nclus, int numClusBefore)
{
	// In this function we know in advance that there is 2 objects at least
	lstFusptr lt = NULL, aux_lt = NULL, lt1=NULL, lt2=NULL;
	fusObptr aux = NULL;
	int i, flag=0, borrar, auxnumC ;

	// First Stage: Extraction of each group and linked list creation
	for (i=1; i<=nclus; i++)
	{
		if ((aux = fillStateObject(Xp, i) )!=NULL)
		{
			//std::cout << "Linking object " << i << std::endl;
			aux_lt = (lstFusptr)calloc(1,sizeof(lstFusptr));
			if(lt==NULL)
			{
				aux_lt->ob = aux;
				aux_lt->prev = NULL;
				aux_lt->next = NULL;
				lt = aux_lt;
			}
			else
			{
				lt->next = aux_lt;   //creo que esta linea esta de mas
				aux_lt->prev = lt;
				lt = lt->next;
				lt->ob = aux;
				lt->next = NULL;
			}
			if (i>numClusBefore)
				flag = nclus;  // New objects will be added
		}
		else
		{
			std::cout << "Problems for memory allocation of object, ANY FUSION REALIZED " << i;
			return 0;
		}
	}
	// Second Stage: Test of proximity and similar velocities
	lt = findHeadFusionList(lt);
	lt1 = lt;
	lt2 = lt->next;
	//  std::cout << " nclus " << nclus << std::endl;
	while(lt1->next != NULL)
	{
		if((fusionTwoObjects(lt1->ob, lt2->ob)) != 0 )
		{
			nclus--;
			//  std::cout << " dentro del lazo, nclus " << nclus << std::endl;
			flag = nclus; // There is a new different object
			// Look for the label of new object (smaller index)
			if (lt1->ob->numClus < lt2->ob->numClus){
				auxnumC = lt1->ob->numClus;
				borrar =lt2->ob->numClus;
			}
			else{
				auxnumC = lt2->ob->numClus;
				borrar = lt1->ob->numClus;
			}
			// Change the Xp jmath matrix
			renameMatrixObjectIndex (Xp, borrar, auxnumC);
			//   std::cout << " Finish rename matrix " << std::endl;
			// Verify if there are more than 2 objects in the list
			//    if (lt->next->next == NULL || nclus <= 2) Si dejo esta linea no prueba 2 grupos finales
			if (lt->next->next == NULL || nclus < 2)
			{ // Termino
				//	   std::cout<< "Just two groups, no more merging test" << std::endl;
				deleteListObjects(lt);
				break;
			}
			else // More than 2 objects
			{
				//   std::cout << "More than 2 objects " << std::endl;
				//Unlink the second element
				aux_lt = lt2->prev;
				aux_lt->next = lt2->next;
				if (lt2->next != NULL)
				{
					aux_lt = lt2->next;
					aux_lt->prev = lt2->prev;
				}
				//  std::cout << "Adios al primero " << std::endl;
				//now the first one
				if (lt1->next == NULL)
				{
					aux_lt =  lt1->prev;
					aux_lt->next = NULL;
				}
				else {
					aux_lt = lt1->next;
					aux_lt->prev = lt1->prev;
				}
				if (lt1->prev != NULL)
				{
					aux_lt = lt1->prev;
					aux_lt->next = lt1->next;
				}
				else
					lt = aux_lt; //new head
				// std::cout << "Adios al segundo " << std::endl;
				// Free the two merging objects
				freeDoubleMatrix(lt1->ob->v);
				freeDoubleMatrix(lt2->ob->v);
				//  std::cout << "Adios a la liberacion de matrices " << std::endl;
				free(lt1->ob);
				free(lt2->ob);
				// std::cout << "Adios a los lt1 y 2 de ob " << std::endl;
			}
			//Link in the Tail of the list the new element
			if ((aux = fillStateObject(Xp, auxnumC) )!=NULL)
			{
				free(lt1);
				free(lt2);
				//	   std::cout << "Linking in the tail New object " << aux->v[(aux->nP)-1][4] << std::endl;
				//	   std::cout << " Para el nuevo objeto Vel1= " << aux->meanVel[0] << " teta1= " << aux->meanVel[1] << std::endl;
				//	   for (int j=0; j<aux->nP; j++)
				//    std::cout<< " 0= " <<aux->v[j][0] << " 1= " <<aux->v[j][1] << " 2= " <<aux->v[j][2] << " 3= " <<aux->v[j][3] <<" 4= " <<aux->v[j][4] << std::endl;
				aux_lt = (lstFusptr)calloc(1,sizeof(lstFusptr));
				aux_lt->next = NULL;
				aux_lt->ob = aux;
				lt1=lt;
				while(lt1->next != NULL) //busca la cola en la lista
					lt1=lt1->next;
				// Enlaza el nuevo cluster en la cola
				lt1->next = aux_lt;
				aux_lt->prev = lt1;
				// set the new objects for fusion iteration
				lt1 = lt;
				lt2 = lt->next;
			} //calcula NFA
			else
			{
				std::cout << "Not Memory reserved for New fusioned object " << std::endl;
				return 0;
			}
			// Values to lt1 and lt2 for the next iteration
			if(lt->next == NULL)
			{  // There is not more groups
				deleteListObjects(lt);
				//  free(lt1);/home/luzdora/MATHWORKS_R2008B_MAC/bin
				// free(lt2);
				break;
			}
			else if(lt1->prev != NULL){ // There is not the head
				aux_lt = lt1;
				lt1 = lt1->prev;
				if (lt1->next != NULL)  // There is not the tail
				{
					lt1=lt1->next;
					free(aux_lt);
					free(lt2);
					lt2 = lt1->next;
				}
			}
			else
			{
				lt1 = lt;
				lt2 = lt->next;
			}
		}
		else  //Not fusion, just set next value for the iteration
			if (lt2->next == NULL)
			{
				lt1 = lt1->next;
				lt2 = lt1->next;
				/*	 if (lt2->next == NULL)
		   {
		     std::cout << "Adios desde parte no hay fusion solo avanza " << std::endl;
		     deleteListObjects(lt);
		     break;
		     } */
			}
			else
				lt2 = lt2->next;
	} // End initial while

	//  if (flag == 0)
	//	 flag = numClusBefore;
	//  std::cout << " Flag = " << flag << std::endl;
	return flag; // Final number of objects in the linked objects list
}

} //end namespace tracker

/*
int main ()
{
	std::cout << "Este es el main de tracker" << std::endl;


	return 0;
}
 */
