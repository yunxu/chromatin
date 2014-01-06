/*
 * =====================================================================================
 *
 *       Filename:  gen.uniform.samples.cpp
 *
 *    Description:  Generate uniform sample points on sphere
 *
 *        Version:  1.0
 *        Created:  02/12/2013 12:52:52
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yun Xu
 *   Organization:  UIC BIOE
 *  
 * =====================================================================================
 */

#include <iostream>
#include "XYSO3Sequence.h"
#include "XYVector.h"
#include "XYFile.h"
#include <vector>
#include	<stdlib.h>

/* vim settings */
/* :set makeprg=make\ -j\ -C\ ../build2 */
/* \rme ../build2/bin/genuniformsamples  */

int main ( int argc, char *argv[] ) {
  int iNumPoints = 160;
  float fSeglen = 1500;

  CXYSO3Sequence sO3sequence(iNumPoints);
  sO3sequence.SetSO3Sequence();
  CXYMatrix<float>* pMSamplesOrg = new CXYMatrix<float>;
  (*pMSamplesOrg) = sO3sequence.GetSO3Sequence();
  cout << "haha" << iNumPoints <<endl;

  CXYMatrix<float> kMSamplePoints = *pMSamplesOrg;
  kMSamplePoints *= fSeglen;

  CXYFile kFile;
  kFile.WriteMatrix("/tmp/UniformSample.pts","w", kMSamplePoints); 

  cout << int(80/30) << endl;
  delete pMSamplesOrg;
  return EXIT_SUCCESS;
}
