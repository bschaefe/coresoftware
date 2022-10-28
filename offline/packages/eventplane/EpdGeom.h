// Tell emacs that this is a C++ source
//  -*- C++ -*-.


/*************************************
 *
 * \author Mike Lisa
 * \date 4 Jan 2018
 *
 * \description:
 *  class defining geometrical aspects of an EPD Tile
 *  (position of center, RandomPointOnTile(), etc.)
 *
 * The user may pass the PP/TT/SN _or_ the uniqueID to
 *   most functions.  No option to pass EpdHit object,
 *   because we want to avoid StObject-dependence.
 *   Instead, simply use EpdHit::id(),
 *   e.g. RandomPointOnTile(hit->id())
 *
 *
 * Adapted for use in sPHENIX (from STAR)
 *   by Brennan Schaefer
 *   11 July 2022
 *
 * Abbreviations
 *   PP = position;
 *   TT = tilenumber;
 *   SN = southnorth;
 *   SS = 30 deg, super sector of two sectors
 *
 *************************************/

#ifndef EPDGEOM_H
#define EPDGEOM_H

#include <gsl/gsl_const.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <TVector3.h>
#include <string>
#include <fun4all/SubsysReco.h>

class PHCompositeNode;

class EpdGeom// : public SubsysReco
{
 public:

  EpdGeom(const std::string &name = "EpdGeom");

//~EpdGeom() override;
  ~EpdGeom();// override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
//  int Init(PHCompositeNode *topNode) override;
  int Init(PHCompositeNode *);// override;
//  int Init(PHCompositeNode *topNode);

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
//int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */


  double GetZwheel();

  void GetCorners(short uniqueID, int* nCorners, double* x, double* y);

  void GetCorners(short position, short tilenumber, short southnorth, int* nCorners, double* x, double* y);

  void GetCorners(int* nCorners, double* xc, double* yc);

  bool IsInTile(short uniqueID, double x, double y);

  bool IsInTile(short position, short tilenumber, short southnorth, double x, double y);

  bool IsInTile(double x, double y);

  TVector3 RandomPointOnTile(short uniqueID);

  TVector3 RandomPointOnTile(short PP, short TT, short SN);

  TVector3 RandomPointOnTile();

  short Row(short uniqueID);

  short Row(short PP, short TT, short SN);

  short Row();

  void SetPpTtSn(short uniqueID);

  TVector3 TileCenter(short uniqueID);

  TVector3 TileCenter(short PP, short TT, short SN);

  TVector3 TileCenter();


 private:

  short mPP;  /// sector position [1,12]
  short mTT;  /// tile number on sector [1,31]
  short mSN;  /// South/North = -1/+1

  double mPhiCenter[12][31][2];  // PP,TT,SN
  double mRmin[16];              // row
  double mRmax[16];              // row
  double mRave[16];              // row

  gsl_rng* m_RandomGenerator = nullptr;

  unsigned int m_Seed = 0;

  gsl_rng* RandomGenerator() const { return m_RandomGenerator; }


};

#endif // EPDGEOM_H
