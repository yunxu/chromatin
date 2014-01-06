/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
/*
 *  XYSO3Sequence.h
 *  CHROMATIN_SIS_COARSE
 *
 *  Created by Yun Xu on 10/11/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef XYSO3SEQUENCE_H
#define XYSO3SEQUENCE_H

#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include "XYMatrix.h"
#include "XYFile.h"
//#include<stdbool.h>

typedef struct Quaternion {
	double w;
	double x;
	double y;
	double z;
} Quaternion;

typedef struct EulerAngle{
	double theta; // pitch
	double psi;		// roll
	double phi;		// yaw
} EulerAngle;

using namespace std;
class CXYSO3Sequence {
public:
	CXYSO3Sequence();
	CXYSO3Sequence(int iNumPoints);
	~CXYSO3Sequence(void);
	void pix2ang_nest(long nside, long ipix, double *theta, double *phi);
	void mk_pix2xy(int *pix2x, int *pix2y);
	vector<double> find_point(int base_grid, long int point,long int level,long int healpix_point);
	bool hopf2quat(vector < vector <double> > Points);
	int SetSO3Sequence();
	CXYMatrix<float>& GetSO3Sequence();
	void	SetSequenceBase();
	vector< int >& GetSequenceBase();
	int GetNumPoints();
	void quat2vector(Quaternion quat, double oldpoint_xyz[], double newpoint_xyz[]);
	void eulerangle2vector(EulerAngle ea, double oldpoint_xyz[], double newpoint_xyz[]);
	void ang2vec(double theta, double phi, double *vec);
	void printSO3Sequence();
private:
	int m_iNumPoints;
	vector < int > m_vSequence_base;
	vector< vector<double > > m_vPoints;
	CXYMatrix<float> *m_pMSpherepoints;
};


#endif
