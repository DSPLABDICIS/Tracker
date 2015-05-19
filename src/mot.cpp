/*
 * mot.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: luzdora
 */
#include "mot.hpp"

namespace tracker {

using namespace boost_incl;
using namespace klt;

inline double sqr(double x) { return x*x; }


/****************************************************************************************

         verify working dims & pos

 ****************************************************************************************/
void verifyWorkingDims(int workingDims[2], int w, int h)
{
	// set to a square because the object can rotate, which could
	// change the dims of the working area (so it should avoid that)
	if (workingDims[0] > workingDims[1])
		workingDims[1] = workingDims[0];
	else
		workingDims[0] = workingDims[1];

	// verify that the box is smaller than the image !
	/*	if (workingDims[0] > h || workingDims[0] > w)
	if (h < w)
	workingDims[0] = workingDims[1] = h; else
	workingDims[0] = workingDims[1] = w;
	 */
}


void verifyWorkingPos(int workingDims[2], int workingPos[2], int w, int h)
{
	// if not possible to put the object at the center, keep the
	// same dimensions but move the working box
	if (workingPos[0] < 0) workingPos[0] = 0;
	if (workingPos[1] < 0) workingPos[1] = 0;
	if (workingPos[0] + workingDims[0] > w)
		workingPos[0] = w-workingDims[0]; // ok because dims[0] < width
	if (workingPos[1] + workingDims[1] > h)
		workingPos[1] = h-workingDims[1]; // ok because dims[1] < height

}

int depthValidationOfGroups(mat & Xp,
		char *filename,
		int nClusters,
		int pActual)
{
	int i, j, accum, nC_aux, index, aux, elem=0, j1=1, shift=0;
	int promDep;
	int depMin=0, depMax=0, depth;
	nC_aux = nClusters;

	if (nClusters > 0){
		IplImage *img;
		img = cvLoadImage(filename, CV_LOAD_IMAGE_GRAYSCALE);

		for (j = 1; j<=nClusters; j++, j1++)
		{
			elem=0;
			accum=0;
			for (i = (int)Xp.size1()-pActual ; i <(int)Xp.size1(); i ++)
				if (Xp(i,4)==j1){
					// Inicia Validacion de profundidad
					elem++;
					//	index = Xp(i,0)*img->widthStep + Xp(i,1);
					index = Xp(i,1)*img->widthStep + Xp(i,0);
					aux = (unsigned char) img->imageData[index];
					accum += aux;
					// Buscar minima y maxima profundidad
					if (elem==1)
						depMin = depMax = aux;
					else
						(depMin>aux)? depMin = aux : ((depMax < aux)? depMax=aux : depMax = depMax );

					std::cout << " Grupo " << j1 << " Pmin "<< depMin << " PMax "<< depMax << " punto X " << Xp(i,0)<< " punto Y " << Xp(i,1)  << " accum "<< accum << " punto " << i << std::endl;
				}

			depth = ( depMax - depMin) ;

			if (elem != 0)
			{
				promDep = accum/elem;
				//	if ( !((depth >= promDep-10) && (depth <= promDep+10) && (promDep>=0)) )
				if (depth > 25)
				{
					// Cambiar el subindice por que el cluster no tiene los elementos a la misma profundidad
					std::cout << "Grupo fuera de la profundidad " << j << " promDep = " << promDep << std::endl;
					for (i=0; i < (int)Xp.size1(); i ++)
						if (Xp(i,4)==j)
							Xp(i,4)=0;
						else
							Xp(i,4)>j? Xp(i,4)-=1:Xp(i,4)=Xp(i,4);
					nC_aux-=1;
					shift++;
					j1= j-shift;
				}
			}
		}
		std::cout << "Numero de clusters despues de la profundidad " << nC_aux << std::endl;
		cvReleaseImage(&img);
	}
	return nC_aux;
}




/****************************************************************************************
	Calcule Convex boundaries
 ****************************************************************************************/
int TrackedObject::calcConvex()
{
	/** Opencv data transform **/
	int i, hullcount ;

	CvPoint* points = (CvPoint*)malloc( (int)fl.size1() * sizeof(points[0]));
	int* hull = (int*)malloc( (int)fl.size1() * sizeof(hull[0]));
	CvMat point_mat = cvMat( 1,(int)fl.size1() , CV_32SC2, points );
	CvMat hull_mat = cvMat( 1, (int)fl.size1(), CV_32SC1, hull );

	for( i = 0; i < (int)fl.size1(); i++)
	{
		points[i].x = (int)(fl(i,0)+0.5);
		points[i].y = (int)(fl(i,1)+0.5);
		//	  	  std::cout << "Punto num= " <<  i << " x " << points[i].x << " y " <<  points[i].y << std::endl;
	}
	cvConvexHull2( &point_mat, &hull_mat, CV_CLOCKWISE, 0 );
	hullcount = hull_mat.cols;
	conv.resize(hullcount, 2, 0);
	//    std::cout << " Puntos en el objeto " << (int)fl.size1() << " puntos en HULL " <<  hullcount << std::endl;

	for( i = 0; i < hullcount; i++ )
	{
		conv(i,0)= points[hull[i]].x;
		conv(i,1)= points[hull[i]].y;
		//	 	  std::cout << " Cont X " << conv(i,0) << " ContY " << conv(i,1) << std::endl;
	}
	//    std::cout << "END OBJECT HULL, points =  "<< hullcount << std::endl;
	free(points);
	free(hull);
	return hullcount;
}


/****************************************************************************************
		calcState
          This function returns the number of points which are tracked in the object set

 ****************************************************************************************/
int TrackedObject::calcState()
{

	int xpos, ypos, i, numP=0;
	// initializations (bounding box to the opposite of max dims, and object pos
	boundingBox[0] = 1000000; boundingBox[1] = 1000000;
	boundingBox[2] = -1;			boundingBox[3] = -1;
	pos[0] = 0; pos[1] = 0;
	vel[0] = 0; vel[1] = 0;

	// adjust bounding box and calc barycenter
	for( i = 0; i < (int)fl.size1(); i++)
	{
		if (fl(i,4)>=0){
			// calc barycenter
			numP++;
			xpos = (int)(fl(i,0)+0.5);
			ypos = (int)(fl(i,1)+0.5);
			pos[0] += xpos;
			pos[1] += ypos;
			preVel[0] = vel[0];
			preVel[1] = vel[1];
			vel[0] += fl(i,2);
			vel[1] += fl(i,3);


			// try iteratively to enlarge the bounding box  por info (enum rect {x1p=0,y1p,x2p,y2p})
			if (xpos < boundingBox[x1p]) boundingBox[x1p] = xpos;
			if (xpos > boundingBox[x2p]) boundingBox[x2p] = xpos;
			if (ypos < boundingBox[y1p]) boundingBox[y1p] = ypos;
			if (ypos > boundingBox[y2p]) boundingBox[y2p] = ypos;
		}
	}

	// finish to calc the barycenter
	/* SI REGIONES
      if (numP > 4 && flagRegion) // Puntos suficientes para contorno y bandera activada
	{
	  pos[0] =  boundingBox[x1p]+ (boundingBox[x2p] - boundingBox[x1p])/2;
	  pos[1] =  boundingBox[y1p]+ (boundingBox[y2p] - boundingBox[y1p])/2;

	}
      else
	 */	if(numP !=0 )
	 {
		 pos[0] /= numP;
		 pos[1] /= numP;
		 vel[0] /= numP;
		 vel[1] /= numP;
	 }
	 else {
		 std::cout << "There is not points" << std::endl;
		 return 0;
	 }
	 //   std::cout << " centre x " << pos[0] << " centre y " << pos[1] << " velx " << vel[0] << " vely " << vel[1] << std::endl;
	 // calc scope
	 double maxDist = 0;
	 for(int i = 0; i < (int)fl.size1(); i++)
	 {
		 if (fl(i,4) >= 0) {
			 double dist = sqr(fl(i,0)-pos[0]) + sqr(fl(i,1)-pos[1]);

			 if (dist > maxDist)
				 maxDist = dist;
		 }
	 }

	 if (!status) // if not the first time, save scope as the previous Scope
		 prevScope = scope;

	 scope = (int)(2*sqrt(maxDist));  // Then calculate the new scope
	 if (status)  // if first time initialise prevScope with scope value
		 prevScope = scope;

	 /*   if (prevScope > 5)
	{
	  if (scope > prevScope*1.1) scope = (int)(prevScope*1.1);
	  if (scope < prevScope/1.1) scope = (int)(prevScope/1.1);
	}
	if (scope < 5) scope = 5;  */
	 //    std::cout << "scope is " << scope << " PrevScope " << prevScope << std::endl;

	 return numP;
}

/****************************************************************************************

		WriteBordersToImg

 ****************************************************************************************/
void TrackedObject::writeKalmanToImg(cv::Mat &jRGBimg)
{
	//using namespace klt;

	// JFR_PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	// JFR_PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	// JFR_PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");

	// CONVEX BORDERS Si quiero convex border poner condicion con signo ">"



	if ((int)conv.size1()>0 && flagRegion){
		for (int i=0;i<(int)fl.size1() ; i++) {
			WritePointToImg(jRGBimg, (int)(fl(i,0) + 0.5), (int)(fl(i,1) + 0.5), 0,0,255); //Blue
			// Solo para que la region no sea considerada como objeto
			//  fl(i,0)= 0.0;  fl(i,1)= 0.0; fl(i,2)= 0.0; fl(i,3)= 0.0;
		}
		int ix = conv((int)conv.size1()-1,0) , iy = conv((int)conv.size1()-1,1);
		for (int i=0;i<(int)conv.size1() ; i++){
			WriteLineToImg(jRGBimg, ix, iy, conv(i,0),conv(i,1),255,255,0);      // 20,230,220); //CYAN
			WriteLineToImg(jRGBimg, ix, iy+1, conv(i,0),conv(i,1)+1,255,255,0);      // 20,230,220); //CYAN
			WriteLineToImg(jRGBimg, ix, iy-1, conv(i,0),conv(i,1)-1,255,255,0);      // 20,230,220); //CYAN

			ix = conv(i,0);
			iy =conv(i,1);
		}
		for (int i=0;i<(int)conv.size1() ; i++)
			WritePointToImg(jRGBimg, (int)(conv(i,0)), (int)(conv(i,1)), 210,50,205); //fiusha
	}
	// ESTO DIBUJA LOS OBJETOS DETECTADOS SIN CONTORNOOOO!!! QUITAR COMENTARIO PORQUE SINO NO HABRA RESULTADO SI NO ESTA ACTIVADA LA BANDERA SNAKE
	else {

//		for (int i=0;i<(int)fl.size1() ; i++)
//			WritePointToImg(jRGBimg, (int)(fl(i,0) + 0.5), (int)(fl(i,1) + 0.5), 0,0,255); //Blue

		//WritePointToImg(jRGBimg, (int)(pos[0] + 0.5), (int)(pos[1] + 0.5), 255,255,0); //Barycento Pasado

		WriteRectToImg(jRGBimg, boundingBox[x1p]-2, boundingBox[y1p]-2,
				boundingBox[x2p]+2, boundingBox[y2p]+2, 255,255,0);
		WriteRectToImg(jRGBimg, boundingBox[x1p]-3, boundingBox[y1p]-3,
				boundingBox[x2p]+3, boundingBox[y2p]+3, 255,255,0);

		//WritePointToImg(jRGBimg, (int)(estimatedPos[0] + 0.5), (int)(estimatedPos[1] + 0.5), 255,0,255); //Barycento Predicho
	}



	// estimation vert
	/*   WritePointToImg(jRGBimg, estimatedPos[0], estimatedPos[1], 255,255,255);

      //		       Kalman en verde previo y en blanco prediction

	WriteRectToImg(jRGBimg, prevWorkingPos[0],prevWorkingPos[1],
	prevWorkingPos[0]+workingDims[0],
	prevWorkingPos[1]+workingDims[1], 0,255,0);  //verde

	WriteRectToImg(jRGBimg, prevWorkingPos[0]-1,prevWorkingPos[1]-1,
	prevWorkingPos[0]+workingDims[0]+1,
	prevWorkingPos[1]+workingDims[1]+1, 0,255,0);  //verde


	WriteRectToImg(jRGBimg, workingPos[0],workingPos[1],
	workingPos[0]+workingDims[0],
	workingPos[1]+workingDims[1], 255,255,255);  //blanco

	 */
}


/****************************************************************************************
		WriteBordersToImg
 ****************************************************************************************/
void TrackedObject::writeBordersToImg(cv::Mat &jRGBimg, KLT_TrackingContext tc)
{
	// using namespace klt;
	// JFR_PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	// JFR_PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	// JFR_PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");


	// AQUI SOLO METER LO QUE IMPRIME BORDES RECTANGULARES Y PUNTOS
//	for (int i=0;i<(int)fl.size1() ; i++)
//		WritePointToImg(jRGBimg, (int)(fl(i,0) + 0.5), (int)(fl(i,1) + 0.5), 0,0,255); //Blue

	WriteRectToImg(jRGBimg, boundingBox[x1p]-2, boundingBox[y1p]-2,
			boundingBox[x2p]+2, boundingBox[y2p]+2, 255,255,0);
	WriteRectToImg(jRGBimg, boundingBox[x1p]-3, boundingBox[y1p]-3,
			boundingBox[x2p]+3, boundingBox[y2p]+3, 255,255,0);

	// barycenter
	//    WritePointToImg(jRGBimg, pos[0], pos[1], 255,0,255);

}


/****************************************************************************************
			Draw contour in image file
 ****************************************************************************************/
void TrackedObject::writeContourToImg(cv::Mat &jRGBimg, unsigned char *featm)
{
	// using namespace klt;
	// int i, j;

	// JFR_PRECOND(jRGBimg.colorSpace()==JfrImage_CS_BGR || jRGBimg.colorSpace()==JfrImage_CS_RGB, "jRGBimg is not in the RGB color space !");
	// JFR_PRECOND(jRGBimg.depth()==IPL_DEPTH_8U, "jRGBimg has not u8 depth !");
	// JFR_PRECOND(jRGBimg.channels()==3, "jRGBimg has not 3 channels !");

	// DISCRET CONTOUR
	//   for ( i=0;i<(int)contour.size1() ; i++)
	//	WritePointToImg(jRGBimg, contour(i,0), contour(i,1), 255,255,0);
	for (int i=0;i<(int)fl.size1() ; i++)
		WritePointToImg(jRGBimg, (int)(fl(i,0) + 0.5), (int)(fl(i,1) + 0.5), 0,0,255); //Blue

	// CONVEX BORDERS  Si quiero convex border poner condicion con signo ">"
	if ((int)conv.size1()>0){
		int ix = conv((int)conv.size1()-1,0) , iy = conv((int)conv.size1()-1,1);
		for (int i=0;i<(int)conv.size1() ; i++){
			WriteLineToImg(jRGBimg, ix, iy, conv(i,0),conv(i,1),255,255,0);      // 20,230,220); //CYAN
			WriteLineToImg(jRGBimg, ix, iy+1, conv(i,0),conv(i,1)+1,255,255,0);      // 20,230,220); //CYAN
			WriteLineToImg(jRGBimg, ix, iy-1, conv(i,0),conv(i,1)-1,255,255,0);      // 20,230,220); //CYAN

			ix = conv(i,0);
			iy = conv(i,1);
		}
		// Convex Points
		for (int i=0;i<(int)conv.size1() ; i++)
			WritePointToImg(jRGBimg, (int)(conv(i,0)), (int)(conv(i,1)), 210,50,205); //fusha

	}


	/*
      for (i=0; i <jRGBimg.width(); i++ )
	  for (j=0; j < jRGBimg.height(); j++)
	     if (featm[j*jRGBimg.width()+i] == 255)
	         WritePointToImg(jRGBimg, i, j, 255,255,255);
	 */
}



/****************************************************************************************
	                        detectLandscape.
     This function returns:
                      0 :  if at least one point is still in the object set
                     -1 : if there is no more points in the object set.
 ****************************************************************************************/

int TrackedObject::detectLandscape(KLT_FeatureList feat)
{
	int i, trackedCount = 0; // count of non lost features
	double mean=0, stdev=0; // mean and standard deviation of displacement norms
	double meanX=0, stdevX=0; // of x displacements
	double meanY=0, stdevY=0; // of y displacements
	double meanDist=0, stdevDist=0; // of distance from features to barycenter
	double auxX, auxY, xpos, ypos;

	pos[0]=0.0;
	pos[1]=0.0;

	// To count only tracked and in the same direction features
	for( i = 0; i < (int)feat->nFeatures; i++)
	{
		if (feat->feature[i]->val >= 0) // if feature non lost
		{
			// EVITARR esta prueba para las sequencias estaticas

			auxX = (feat->feature[i]->x-fl(i,0))/dt;
			auxY = (feat->feature[i]->y-fl(i,1))/dt;

			if ( (auxX*vel[0]>=0 || (fabs(auxX) > LIM )) && (auxY*vel[1]>=0 || (fabs(auxY) > LIM )))   // the same direction and not so slowly Todo: What?????
				//	if ( (auxX*vel[0]>=0)  || (fabs(auxX) < LIM )  )   // the same direction and not so slowly
			{
				trackedCount++;
				xpos =  feat->feature[i]->x ;
				ypos =  feat->feature[i]->y;
				pos[0] += xpos;
				pos[1] += ypos;
				//std::cout<< "NOOO X= "<< feat->feature[i]->x<< " NOOO Y= " << feat->feature[i]->y<< " Not Landscape..."<<std::endl;

			}
			else {
				feat->feature[i]->val = KLT_LANDSCAPE;
				//std::cout<< "NOOO auxX= " << auxX << " vel[0]= " << vel[0] << " feat(t) " <<feat->feature[i]->x  << " feat(t-1) "<< fl(i,0)  << std::endl;
				//std::cout<< "NOOO auxY= " << auxY << " vel[1]= " << vel[1] << " feat(t) " <<feat->feature[i]->y  << " feat(t-1) "<< fl(i,1)  << std::endl;
			}
		}
	}

	if ( trackedCount != 0)
	{
		pos[0] /= trackedCount;
		pos[1] /= trackedCount;
	}
	else {
		//std::cout << " LANDSCAPE : There is not any feature in this object " << std::endl;
		return -1; // There is not any feature in this object
	}

	// If there are features
	double *displacements = new double[feat->nFeatures]; // norms of displacements
	double *dispX = new double[feat->nFeatures]; // displacements along x
	double *dispY = new double[feat->nFeatures]; // displacements along y
	double *distances = new double[feat->nFeatures]; // distances from features to the barycenter

	for( i = 0; i < feat->nFeatures; i++)
	{
		if (feat->feature[i]->val >= 0) // if feature non lost
		{ // calc displacements
			dispX[i] = feat->feature[i]->x-fl(i,0);
			dispY[i] = feat->feature[i]->y-fl(i,1);
			displacements[i] = sqrt(sqr(dispX[i])+sqr(dispY[i]));
			distances[i] = sqrt(sqr(feat->feature[i]->x-pos[0]) + sqr(feat->feature[i]->y-pos[1]));
			//  std::cout << "feature " << i << " x = " << feat->feature[i]->x << " ; y = " << feat->feature[i]->y << " dispX = " << dispX[i] << " dispY = " << dispY[i] << " disp = " << displacements[i] << " dist = " << distances[i] << std::endl;
			//  update means
			mean += displacements[i];
			meanX += dispX[i];
			meanY += dispY[i];
			meanDist += distances[i];
		}
	}

	// finish calc means
	mean  /= trackedCount;
	meanX /= trackedCount;
	meanY /= trackedCount;
	meanDist /= trackedCount;

	// calc standard devs
	for( i = 0; i < feat->nFeatures; i++)
	{
		if (feat->feature[i]->val < 0) continue;
		stdev  += sqr(displacements[i]-mean);
		stdevX += sqr(dispX[i]-meanX);
		stdevY += sqr(dispY[i]-meanY);
		stdevDist += sqr(distances[i]-meanDist);
	}

	stdev /= trackedCount; stdev = sqrt(stdev);    // if (stdev < 2) stdev = 2;
	stdevX /= trackedCount; stdevX = sqrt(stdevX); // if (stdevX < 2) stdevX = 2;
	stdevY /= trackedCount; stdevY = sqrt(stdevY); // if (stdevY < 2) stdevY = 2;
	stdevDist /= trackedCount; stdevDist = sqrt(stdevDist); // if (stdevDist < 2) stdevDist = 2;
	//   std::cout << "mean " << mean << ", stdev " << stdev << ", meanX " << meanX << ", stdevX " << stdevX << ", meanY " << meanY << ", stdevY " << stdevY << ", meanDist " << meanDist << ", stdevDist " << stdevDist << std::endl;

	// suppress bad features (in the landscape) and update new
	// barycenter pos (without suppressed features)
	double coeff;
	boundingBox[0] = 1000000; boundingBox[1] = 1000000;
	boundingBox[2] = -1;			boundingBox[3] = -1;
	pos[0] = 0.0; pos[1] = 0.0;
	trackedCount = 0;
	for(i = 0; i < feat->nFeatures; i++)
	{
		if (feat->feature[i]->val < 0) continue; // if non lost

		// if (feat->feature[i]->val > 0) coeff = 1.2; else coeff = 2;
		coeff = 2.6;
		// first test on displacements
		/*	  if ((dispX[i] > (meanX+coeff*stdevX + 0.5)) || (dispX[i] < (meanX-coeff*stdevX-0.5)) ||
	    (dispY[i] > (meanY+coeff*stdevY + 0.5)) || (dispY[i] < (meanY-coeff*stdevY-0.5)))  */
		if ((dispX[i] > (meanX+coeff*stdevX + 0.5)) || (dispX[i] < (meanX-coeff*stdevX-0.5)))
		{
			feat->feature[i]->val = KLT_LANDSCAPE;
			//std::cout << "feat " << i << " in landscape (dispX " << dispX[i] << ", dispY " << dispY[i] << ") LimitX " << meanX+coeff*stdevX+0.5 << " inf " << meanX-coeff*stdevX-0.5 << ") LimitY " << meanY+coeff*stdevY+0.5 << " inf " << meanY-coeff*stdevY-0.5 << std::endl;
			continue;
		}

		// second test on displacement norm : if displacement
		// of one feature is too large or too small related to
		// others, suppress it  (scope/3)


		if ((displacements[i] > (mean+stdev/2)*coeff && displacements[i] > scope) ||
				(displacements[i] < (mean-stdev/2)/coeff))
		{
			feat->feature[i]->val = KLT_LANDSCAPE;
			//std::cout << "feature " << i << " in landscape (displacement " << displacements[i] << ") > que " << (mean+stdev/2)*coeff  << " menor que " <<(mean-stdev/2)/coeff << " scope " << scope << std::endl;
			continue;
		}

		/*
	  // third test on position : if a feature is really far from
			  // the others, suppress it
	  //coeff = 1.5;
	  	  coeff = 2.0;
	  if ((distances[i] > (meanDist+stdevDist/2)*coeff))
	    //  if (distances[i] > (meanDist+stdevDist*coeff))
	    {
	      feat->feature[i]->val = KLT_LANDSCAPE;
	      std::cout << "feature " << i << " in landscape (distance " << distances[i] << ") vs " << (meanDist+stdevDist*coeff) << std::endl;
	      continue;
	    }
		 */

		// update pos
		trackedCount++;
		pos[0] += (int)(feat->feature[i]->x+0.5);
		pos[1] += (int)(feat->feature[i]->y+0.5);
		// try iteratively to enlarge the bounding box
		if ((int)(feat->feature[i]->x+0.5) < boundingBox[x1p]) boundingBox[x1p] = (int)(feat->feature[i]->x+0.5);
		if ((int)(feat->feature[i]->x+0.5) > boundingBox[x2p]) boundingBox[x2p] = (int)(feat->feature[i]->x+0.5);
		if ((int)(feat->feature[i]->y+0.5) < boundingBox[y1p]) boundingBox[y1p] = (int)(feat->feature[i]->y+0.5);
		if ((int)(feat->feature[i]->y+0.5) > boundingBox[y2p]) boundingBox[y2p] = (int)(feat->feature[i]->y+0.5);

	}

	if (trackedCount != 0)
	{
		pos[0] /= trackedCount;
		pos[1] /= trackedCount;
	}

	delete[] displacements;
	return 0;
}

/*

    void TrackedObject::detectLandscape(KLT_FeatureList feat)
    {
      int i, trackedCount = 0; // count of non lost features
      double *displacements = new double[feat->nFeatures]; // norms of displacements
      int *dispX = new int[feat->nFeatures]; // displacements along x
      int *dispY = new int[feat->nFeatures]; // displacements along y
      double *distances = new double[feat->nFeatures]; // distances from features to the barycenter

      double mean=0, stdev=0; // mean and standard deviation of displacement norms
      double meanX=0, stdevX=0; // of x displacements
      double meanY=0, stdevY=0; // of y displacements
      double meanDist=0, stdevDist=0; // of distance from features to barycenter

      int xpos, ypos;

      pos[0]=0;
      pos[1]=0;
      for( i = 0; i < (int)feat->nFeatures; i++)
	{
	  if ( feat->feature[i]->val<0)
	    continue;
	  // calc barycenter
	  trackedCount++;
	  xpos = (int)( feat->feature[i]->x + 0.5);
	  ypos = (int)( feat->feature[i]->y + 0.5);
	  pos[0] += xpos;
	  pos[1] += ypos;
	}

      if ( trackedCount != 0)
	{
	  pos[0] /= trackedCount;
	  pos[1] /= trackedCount;
	}


	for( i = 0; i < feat->nFeatures; i++)
	{
	  if (feat->feature[i]->val < 0) continue; // if feature non lost
	  //  trackedCount++;

	  // calc displacements
	  dispX[i] = (int)(feat->feature[i]->x-fl(i,0));
	  dispY[i] = (int)(feat->feature[i]->y-fl(i,1));
	  displacements[i] = sqrt(sqr(dispX[i])+sqr(dispY[i]));
	  distances[i] = sqrt(sqr(feat->feature[i]->x-pos[0]) + sqr(feat->feature[i]->y-pos[1]));
	  //	  std::cout << "feature " << i << " x = " << feat->feature[i]->x << " ; y = " << feat->feature[i]->y << " dispX = " << dispX[i] << " dispY = " << dispY[i] << " disp = " << displacements[i] << " dist = " << distances[i] << std::endl;

	  // update means
	  mean += displacements[i];
	  meanX += dispX[i];
	  meanY += dispY[i];
	  meanDist += distances[i];
	}

      // finish calc means
      mean /= trackedCount;
      meanX /= trackedCount;
      meanY /= trackedCount;
      meanDist /= trackedCount;

      // calc standard devs
      for( i = 0; i < feat->nFeatures; i++)
	{
	  if (feat->feature[i]->val < 0) continue;
	  stdev += sqr(displacements[i]-mean);
	  stdevX += sqr(dispX[i]-meanX);
	  stdevY += sqr(dispY[i]-meanY);
	  stdevDist += sqr(distances[i]-meanDist);
	}

      stdev /= trackedCount; stdev = sqrt(stdev); if (stdev < 2) stdev = 2;
      stdevX /= trackedCount; stdevX = sqrt(stdevX); if (stdevX < 2) stdevX = 2;
      stdevY /= trackedCount; stdevY = sqrt(stdevY); if (stdevY < 2) stdevY = 2;
      stdevDist /= trackedCount; stdevDist = sqrt(stdevDist); if (stdevDist < 2) stdevDist = 2;
      //   std::cout << "mean " << mean << ", stdev " << stdev << ", meanX " << meanX << ", stdevX " << stdevX << ", meanY " << meanY << ", stdevY " << stdevY << ", meanDist " << meanDist << ", stdevDist " << stdevDist << std::endl;

      // suppress bad features (in the landscape) and update new
      // barycenter pos (without suppressed features)
      double coeff;
      boundingBox[0] = 1000000; boundingBox[1] = 1000000;
      boundingBox[2] = -1;			boundingBox[3] = -1;
      pos[0] = 0; pos[1] = 0;
      trackedCount = 0;
      for(i = 0; i < feat->nFeatures; i++)
	{
	  if (feat->feature[i]->val < 0) continue; // if non lost

	  if (feat->feature[i]->val > 0) coeff = 1.2; else coeff = 2;

	  // first test on displacements
	  if ((dispX[i] > (meanX+coeff*stdevX)) || (dispX[i] < (meanX-coeff*stdevX)) ||
	      (dispY[i] > (meanY+coeff*stdevY)) || (dispY[i] < (meanY-coeff*stdevY)))
	    {
	      feat->feature[i]->val = KLT_LANDSCAPE;
	      std::cout << "feature " << i << " in landscape (dispX " << dispX[i] << ", dispY " << dispY[i] << ")" << std::endl;
	      continue;
	    }

	  // second test on displacement norm : if displacement
	  // of one feature is too large or too small related to
	  // others, suppress it


	//  if ((displacements[i] > (mean+stdev/2)*coeff && displacements[i] > scope/3) ||
	 //     (displacements[i] < (mean-stdev/2)/coeff))
	 //   {
	 //     feat->feature[i]->val = KLT_LANDSCAPE;
	 //     std::cout << "feature " << i << " in landscape (displacement " << displacements[i] << ")" << std::endl;
	 //     continue;
	 //   }


	  // third test on position : if a feature is really far from
			  // the others, suppress it
	  //coeff = 1.5;
	  coeff = 2.0;
	  if ((distances[i] > (meanDist+stdevDist/2)*coeff))
	    //  if (distances[i] > (meanDist+stdevDist*coeff))
	    {
	      feat->feature[i]->val = KLT_LANDSCAPE;
	      std::cout << "feature " << i << " in landscape (distance " << distances[i] << ") vs " << (meanDist+stdevDist*coeff) << std::endl;
	      continue;
	    }


	  // update pos
	  trackedCount++;
	  pos[0] += (int)(feat->feature[i]->x+0.5);
	  pos[1] += (int)(feat->feature[i]->y+0.5);
	  // try iteratively to enlarge the bounding box
	  if ((int)(feat->feature[i]->x+0.5) < boundingBox[x1p]) boundingBox[x1p] = (int)(feat->feature[i]->x+0.5);
	  if ((int)(feat->feature[i]->x+0.5) > boundingBox[x2p]) boundingBox[x2p] = (int)(feat->feature[i]->x+0.5);
	  if ((int)(feat->feature[i]->y+0.5) < boundingBox[y1p]) boundingBox[y1p] = (int)(feat->feature[i]->y+0.5);
	  if ((int)(feat->feature[i]->y+0.5) > boundingBox[y2p]) boundingBox[y2p] = (int)(feat->feature[i]->y+0.5);

	}

      if (trackedCount != 0)
	{
	  pos[0] /= trackedCount;
	  pos[1] /= trackedCount;
	}

      delete[] displacements;
    }

 */

} // end tracker

