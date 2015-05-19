/*
 * mergingObjs.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: luzdora
 */
#include "tracker.hpp"


  namespace tracker {

  using namespace boost_incl;
  using namespace klt;
 // using namespace filter;

    void initialiseStructObject(mat obj,  objMptr ob1, int numC, int nP, int ptsLastIm )
    {
      int i, l, j, dim;
      ob1->nP = nP;

      //ob1->numClus = numC;

      dim = obj.size2();
      ob1->o->fl.resize(nP,dim,0);

      for( i = (int)obj.size1() - ptsLastIm, l=0; i < (int)obj.size1(); i++)
	if (obj(i,4) == numC)
	  {
	    for (j =0; j < dim; j++)
	      ob1->o->fl(l,j) = obj(i,j);
	      //std::cout << " Objeto " << ob1->o->fl(l,4) << " Gradiente " << ob1->o->fl(l,5) << std::endl;
	      //std::cout << " Posicion " << ob1->o->fl(l,6) << "\t" << ob1->o->fl(l,7)<< ob1->o->fl(l,6) << "\t" << ob1->o->fl(l,8) << std::endl;
	    ob1->o->fl(l,4) = 0; // Feature tracked
	    l++;
	  }
    }

    // This function needs a previous verification of ptsLastIm >0

    int countNumberOfObjectPoints(mat _Xp,int ptsLastIm, int nc)
    {
      int i, l;

      for( i = _Xp.size1() - ptsLastIm, l=0; i < (int)_Xp.size1(); i++)
	{
	  if (_Xp(i,4) == nc) // Compare if the point belongs to cluster nc
	    l++;
	  //std::cout<< "punto  " << i << " cluster " << _Xp(i,4) << std::endl;
	}
      //     std::cout<< "Cluster:  " << nc << " Number of points " << l << std::endl;
      return l;
    }

    void renameMatrixObjectIndex (mat &Xp, int name, int rename)
    {
      int i;
      //  std::cout << "Cambiar name = " << name << "Por = " << rename << std::endl;

      for (i=0; i<(int)Xp.size1(); i++)
	if (Xp(i,4)==name){
	  Xp(i,4) = rename;
	}
    }

    fusObptr fillStateObject (mat & obj, int nC)
    {
      fusObptr ob1;

      ob1 = (fusObptr)calloc(1,sizeof(fusionObs));
      if ((ob1->nP = pointsToStructure(obj, ob1, nC)) == 0)
	{
	  std::cout << "There is not points for this cluster" << std::endl;
	  free(ob1);
	  return NULL;
	}
      ob1->numClus = nC;
      ob1->state[0] =  maximalPos(ob1->v, ob1->nP, 0);
      ob1->state[1] =  minimalPos(ob1->v, ob1->nP, 0);
      ob1->state[2] =  maximalPos(ob1->v, ob1->nP, 1);
      ob1->state[3] =  minimalPos(ob1->v, ob1->nP, 1);
      ob1->meanVel[0] = meanOfPoints(ob1->v, ob1->nP, 2);
      ob1->meanVel[1] = meanOfPoints(ob1->v, ob1->nP, 3);
      //  ob1->meanVel[1] = meanOfOrientation(ob1->v, ob1->nP);
      return ob1;
    }

    int fusionTwoObjects(fusObptr ob1, fusObptr ob2)
    {
      //  std::cout << " Fusion Test "<< ob1->numClus << " and " << ob2->numClus << std::endl;
      //  std::cout << " Vel1= " << ob1->meanVel[0] << " Vel2= " <<  ob2->meanVel[0] << " teta1= " << ob1->meanVel[1] << " Teta2= " << ob2->meanVel[1] << std::endl;
	 //	 if ( fabs(ob1->meanVel[0] - ob2->meanVel[0]) < 20 )  // 2.0
	   //   if( (fabs(ob1->meanVel[1] - ob2->meanVel[1]) < 12) || (fabs(ob1->meanVel[1] - ob2->meanVel[1]) > 348) )
      if ( fabs(fabs(ob1->meanVel[0]) - fabs(ob2->meanVel[0])) < LIM && ((ob1->meanVel[0]*ob2->meanVel[0])>0) )  // 2.0
	   //  if ( fabs(fabs(ob1->meanVel[1]) - fabs(ob2->meanVel[1])) < 20 && (ob1->meanVel[1]*ob2->meanVel[1])>0 )  // 2.0
	   if ( fabs(fabs(ob1->meanVel[1]) - fabs(ob2->meanVel[1])) < 10 )
	     {
	       //    std::cout << "Similar velocity and orientation" << std::endl;
	       //second Test : Proximity Position
	       double intX, intY;
	       intX = cluster::intersectRegion(ob1->state[0]+10.0, ob1->state[1]-10.0, ob2->state[0]+10.0, ob2->state[1]-10.0);
	       intY = cluster::intersectRegion(ob1->state[2]+10.0, ob1->state[3]-10.0, ob2->state[2]+10.0, ob2->state[3]-10.0);
	       if ((intX > 0) && (intY > 0))
		 {
		   //	   std::cout << "There is an intersectionin positions: FUSION IS POSSIBLE" << std::endl;
		   return 1; // Intersection in regions
		 }
	       //    else
	       //	 std::cout<< "They are not close " << std::endl;
	     }
      //   else
      //     std::cout<< "Not similar orientation" << std::endl;
      //	 else
      //	   std::cout<< "Not similar velocity";

	 return 0;
    }

    void readMascara (bool *mascara, float tamanio)
    {
      FILE *fp;
      // int i;
      char fname[100] = "/tmp/mask.bin";

      fp = fopen(fname , "rb");
      if (fp == NULL)
        printf(" SNAKE: Can't open file '%s' for reading", fname);

      fread(mascara, tamanio, sizeof(bool), fp);
      //  printf("\n Tamanio MASCARA = %d \n ",sizeof(mascara) );
      /*  for (i = 0 ; i < 640*480 ; i++)  {
        printf(" %d ", (int)mascara[i]);
        } */
      fclose(fp);
    }

    double meanOfPoints(double **_Xpol, int nP, int coor)
    {
      int i;
      double aux = 0.0;

       for (i=0; i<nP; i++)
	 aux+=_Xpol[i][coor];

       aux/=(double)nP;

       return aux;
    }


    double meanOfOrientation(double **_Xpol, int nP)
    {
      int i;
      double min, max;
      double aux = 0.0;

      min = minimalPos(_Xpol, nP, 3);
      max = maximalPos(_Xpol, nP, 3);

      if (max-min > 345){
	double r = 360.0-max;
	min += r;
	aux = max + min/2;
	if (aux >360)
	  aux-=360;
      }
      else{
	for (i=0; i<nP; i++)
	  aux+=_Xpol[i][3];

	aux/=(double)nP;
      }

       return aux;
    }


    double minimalPos(double **_Xpol, int nP, int coor)
    {
      int i;
      double minX = _Xpol[0][coor];
      for (i=1; i<nP; i++)
	if (_Xpol[i][coor] < minX)
	  minX = _Xpol[i][coor];
      return minX;
    }

    double maximalPos(double **_Xpol, int nP, int coor)
    {
      int i;
      double maxX = _Xpol[0][coor];
      for (i=1; i<nP; i++)
	if (_Xpol[i][coor] > maxX)
	  maxX = _Xpol[i][coor];
      return maxX;
    }

    void orderNumberOfObjects(mat &_Xp, int nClus)
    {
      int *t;
      int i, j;

      t = (int*) calloc(nClus,sizeof(int));
      t[0]=1;
      for (i=1; i<nClus; i++)
	{
	  t[i] = 1000;
	  for (j=0; j< (int)_Xp.size1(); j++)
	    if((_Xp(j,4)>t[i-1]) && (_Xp(j,4)<t[i]))
	      t[i]=_Xp(j,4);
	}

      for(i=0; i< (int)_Xp.size1(); i++)
	for (j=0; j< nClus; j++)
	  if(_Xp(i,4) == t[j]){
	    _Xp(i,4) = j+1;
	    //	    std::cout << " Num cluster del punto " << i << " = "  << _Xp(i,4) << std::endl;
	    break;
	  }
      free(t);
    }// end function

    int pointsToStructure(mat _Xp, fusObptr &aux, int nc)
    {
      int i, l, count;
      count = countNumberOfObjectPoints(_Xp, (int)_Xp.size1(), nc);
      if (count >0 )
	{
	  aux->v =  cluster::matrixDoubleReserve(count,5);
	  for( i = 0, l=0; i < (int)_Xp.size1(); i++)
	    {
	      if ((int)_Xp(i,4) == nc)
		{
		  aux->v[l][0] = _Xp(i,0);
		  aux->v[l][1] = _Xp(i,1);
		  aux->v[l][2] = _Xp(i,2);
		  aux->v[l][3] = _Xp(i,3);
		  aux->v[l][4] = _Xp(i,4); // tracked
		  l++;
		}
	    }
	  //	  std::cout << "End extraction " << l << " points of the object " << nc << std::endl;
	  return count;
	}
      else
	return 0; //there is not points for that cluster
    }


    lstObPtr findHeadOfList(lstObPtr list)
    {
      while(list->prev!=NULL)
	list=list->prev;
      return list;
    }

    lstFusptr findHeadFusionList(lstFusptr list)
    {
      while(list->prev!=NULL)
	list=list->prev;
      return list;
    }

    lstObPtr findTailOfList(lstObPtr list)
    {
      while(list->next!=NULL)
	list=list->next;
      return list;
    }

   void freeDoubleMatrix(double **m)
    {
      free(m[0]);
      free(m);
    }

    void deleteListObjects (lstFusptr & _eG)
    {
      _eG = findHeadFusionList(_eG);
      //      std::cout << "Encontro cabeza " << std::endl;
      while(_eG != NULL) {
	//	freeDoubleMatrix(_eG->ob->coord);
	freeDoubleMatrix(_eG->ob->v);
	free(_eG->ob);
	if (_eG->next != NULL) {
	  _eG = _eG->next;
	  free(_eG->prev);
	}
	else {
	  free(_eG);
	  _eG = NULL;
	}
      }
    }

    void eraseListObjects (lstObPtr & _eG)
    {
      while(_eG != NULL) {
	delete (*(_eG->ob->o)).kalman;
	delete _eG->ob->o;
	free(_eG->ob);
	if (_eG->next != NULL) {
	  _eG = _eG->next;
	  free(_eG->prev);
	}
	else {
	  free(_eG);
	  _eG = NULL;
	}
      }
      //   std::cout << "Free Object list " << std::endl;
    }

    void shiftNumberOfObject(lstObPtr &list, int eraser)
    {
      while(list->next!=NULL)
	{
	  if (list->ob->numClus > eraser)
	    list->ob->numClus-=1;
	  list=list->next;
	}
    }

    void saveObjectsOfLinkedList(lstObPtr &lt, int nOb)
    {
      int i;
      for (i=1; i<= nOb; i++)
	lt = lt->next;
      lt->prev->next = NULL; // close the linked list
      lt->prev = NULL; //new head of list

      // Erase Objects
      eraseListObjects(lt);
    }

    void addObjectsOfLinkedList(lstObPtr &lt, int nOb)
    {
      lstObPtr  aux_lt= NULL;
      int i;

      lt = findTailOfList(lt);

      for (i =1 ; i<= nOb; i++)
	{
	  aux_lt = (lstObPtr)calloc(1,sizeof(lstObPtr));
	  aux_lt->ob = (objMptr)calloc(1,sizeof(objM));

	  if (aux_lt->ob != NULL)
	    aux_lt->ob->o = new TrackedObject;
	  else
	    std::cout << "Problems, No memory reserved for new linked Object " << i << std::endl;
	  lt->next = aux_lt;   //creo que esta linea esta de mas
	  aux_lt->prev = lt;
	   lt = aux_lt;
	   lt->next = NULL;
	}
    }

    int extractPointsOfObject(mat &_Xp, mat & obj, int nc)
    {
      int i, l, count=0;

      for(i = 0; i < (int)_Xp.size1(); i++)
	{
	  if ((int)_Xp(i,4)==nc)
	    count++;
	}
      //    std::cout << "Points founded : " << count << " to cluster:  " << nc << std::endl;
      if (count >0 )
	{
	  mat aux(count,5);

	  for( i = 0, l=0; i < (int)_Xp.size1(); i++)
	    {
	      if ((int)_Xp(i,4) == nc)
		{
		  aux (l,0) = _Xp(i,0);
		  aux (l,1) = _Xp(i,1);
		  aux (l,2) = _Xp(i,2);
		  aux (l,3) = _Xp(i,3);
		  aux (l,4) = 0; // tracked
		  l++;
		}
	    }
	  obj = aux;  // Illegal asignation for memory free problems ??? smth like that
	  return 0;
	}
      else
	return 1; //there is not points for that cluster
    }
  }



