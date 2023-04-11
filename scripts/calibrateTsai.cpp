/*!
 \file calibrateTsai.cpp
 \brief Camera - robot extrinsic calibration for the St�ubli robots
 *      Extension of the ViSP library example calibrateTsai ($Id:
 *      calibrateTsai.cpp 3619 2012-03-09 17:28:57Z fspindle $)
 *
 * See https://home.mis.u-picardie.fr/~g-caron/fr/index.php?page=7 for more information.
 *
 * This software was developed at:
 * MIS - UPJV
 * 33 rue Saint-Leu
 * 80039 AMIENS CEDEX
 * France
 *
 * ViSP software is free software and is modified and redistributed
 * under the terms of the GNU General Public License
 * ("GPL") version 2 as published by the Free Software Foundation.
 *
 * See http://www.irisa.fr/lagadic/visp/visp.html for more information.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * Description:
 * Tsai calibration example to estimate hand to eye transformation on a St�ubli robot.
 *
 * Authors:
 * Guillaume Caron (update)
 * Fabien Spindler (original author)
 *
 *****************************************************************************/


/*!
  \example calibrateTsai.cpp
  \brief Example of Tsai calibration to estimate extrinsic camera parameters, ie hand-eye homogeneous transformation.

*/

#include <visp/vpDebug.h>
#include <visp/vpParseArgv.h>
#include <visp/vpIoTools.h>
#include <visp/vpCalibration.h>
#include <visp/vpExponentialMap.h>
#include <stdio.h>
#include <sstream>
#include <iomanip>

#include <fstream>

int main()
{
  // We want to calibrate the hand to eye extrinsic camera parameters from 6 couple of poses: cMo and wMe
  const int N = 6; //4 suffisent
  // Input: six couple of poses used as input in the calibration proces
  std::vector<vpHomogeneousMatrix> v_cMo; // eye (camera) to object transformation. The object frame is attached to the calibrartion grid
    v_cMo.reserve(N);
  std::vector<vpHomogeneousMatrix> v_wMe; // world to hand (end-effector) transformation
    v_wMe.reserve(N);
  // Output: Result of the calibration
  vpHomogeneousMatrix eMc, cMe, eMct; // hand (end-effector) to eye (camera) transformation

  std::ifstream fic_cMo("../data/cMo.txt"), fic_wMe("../data/wMe_UR10.txt");

  //chargement des cMo (unit� : m)
  unsigned i = 0;
  while((i < N) && !fic_cMo.eof())
  {
      v_cMo.push_back(vpHomogeneousMatrix());
	  for(unsigned l = 0 ; l < 4 ; l++)
		  for(unsigned c = 0 ; c < 4 ; c++)
			fic_cMo >> v_cMo[i][l][c];
      
      //v_cMo[i] = v_cMo[i].inverse(); //For Nathan Crombez only
      
	  std::cout << "cMo[" << i <<"] = " << v_cMo[i] << std::endl;

	i++;
  }
  fic_cMo.close();

  //chargement des poses Staubli (unit�s : mm et �)
  i = 0;
  vpHomogeneousMatrix eMw;
  vpRzyxVector rzyx;
  vpRotationMatrix eRw, cRe;
  vpTranslationVector etw, cte;
  vpColVector wE(3), eE(3), eC(3);
  double alpha;
  while((i < N) && !fic_wMe.eof())
  {
	for(unsigned l = 0 ; l < 3 ; l++)
	{
		fic_wMe >> wE[l];//T[l][3];
		//T[l][3] *= 0.001;
		wE[l] *= 0.001;
	}

	//std::cout << "p[" << i <<"] = " << p << std::endl;

	for(unsigned l = 0 ; l < 3 ; l++)
	{
		fic_wMe >> alpha;
		rzyx[2-l] = -alpha * M_PI/180.0;
		// rzyx[2-l] = alpha;
	}

	std::cout << "rzyx[" << i <<"] = " << rzyx.t() << std::endl;

	eRw.buildFrom(rzyx);
	eE = eRw*wE;
	for(unsigned l = 0 ; l < 3 ; l++)
		etw[l]=-eE[l];

	eMw.buildFrom(etw, eRw);

      v_wMe.push_back(eMw.inverse());
      // v_wMe.push_back(eMw);

	std::cout << "wMe[" << i <<"] = " << v_wMe[i] << std::endl;

	i++;
  }
  fic_wMe.close();

    
i = 0;
std::cout << "verif : " << std::endl;
for(;i < N ; i++)
{
    eMw = v_wMe[i].inverse();
  eMw.extract(etw);
  eMw.extract(eRw);

  std::cout << "pose Staubli : " << std::endl;
  wE = eRw.t()*(-etw);
  std::cout << 1000.0*wE.t() << " ";
  rzyx.buildFrom(eRw);
  for(unsigned l = 0 ; l < 3 ; l++)
	std::cout << (180.0/M_PI)*rzyx[2-l] << " ";
  std::cout  << std::endl;
}

  // Compute the eMc hand to eye transformation from six poses
  // - cMo[6]: camera to object poses as six homogeneous transformations
  // - wMe[6]: world to hand (end-effector) poses as six homogeneous transformations
  vpCalibration::calibrationTsai(v_cMo, v_wMe, eMc) ;

  std::cout << std::endl << "Output: hand to eye calibration result: eMc estimated " << std::endl ;
  std::cout << eMc << std::endl ;
    
  cMe = eMc.inverse();
   
  cMe.extract(cte);
  cMe.extract(cRe);

  std::cout << "pose Staubli : X(mm), Y(mm), Z(mm), rX(deg), rY(deg), rZ(deg)" << std::endl;
  eC = cRe.t()*(-cte);
  std::cout << 1000.0*eC.t() << " ";
  rzyx.buildFrom(cRe);
  for(i = 0 ; i < 3 ; i++)
	std::cout << -(180.0/M_PI)*rzyx[2-i] << " ";
  std::cout  << std::endl;
  std::cout << "pose Staubli : X(mm), Y(mm), Z(mm), rX(rad), rY(rad), rZ(rad)" << std::endl;
  std::cout << 1000.0*eC.t() << " ";
  for(i = 0 ; i < 3 ; i++)
	std::cout << -rzyx[2-i] << " ";
  std::cout  << std::endl;
  return 0 ;
}