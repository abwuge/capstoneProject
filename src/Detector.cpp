#include "Detector.h"

#include <algorithm>
#include <thread>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Point2D.h>
#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "MemoryPool.h"

Detector::Detector(const std::vector<ScintillatorCounters> &scintillatorCounters, const ROOT::Math::XYZVector &B)
    : dScintillatorCounters(scintillatorCounters), dB(B) {
  std::sort(
      this->dScintillatorCounters.begin(),
      this->dScintillatorCounters.end(),
      [](const ScintillatorCounters &a, const ScintillatorCounters &b) { return a.getLocation() < b.getLocation(); }
  );
}

Detector::~Detector() {}

HitsData *Detector::particleHitData(
    const Particle &particleOrignal,
    const bool      enableEnergyLoss,
    const bool      enableEnergyLossFluctuation,
    const HitsData *measuresData
) const {
  // use MemoryPool to avoid futex_wait in std::mutex
  HitsData *hitsData = MemoryPool::getDetector_particleHitData_hitsData(this->dScintillatorCounters.size());
  Particle *particle = MemoryPool::getDetector_particleHitData_particle();

  *particle = particleOrignal;

  double       time = 0, propagationLength = 0;
  const double initialZ = this->dScintillatorCounters.front().getLocation();
  for (const ScintillatorCounters &scintillatorCounter : this->dScintillatorCounters) {
    constexpr double conversionFactor__c2cm_ns = TMath::Ccgs() * 1e-9; // conversion factor from c to cm/ns

    const ROOT::Math::XYZPoint  &position = particle->getPosition();
    const ROOT::Math::XYZVector &velocity = particle->getVelocity() * conversionFactor__c2cm_ns; // velocity in cm/ns
    const int                    charge   = particle->getCharge();

    double               deltaTime = 0, deltaPropagatedLength = 0;
    ROOT::Math::XYZPoint hitPosition;

    auto linearPropagate = [&]() -> bool {
      if (velocity.Z() <= 0) {
        if (Config::enableWarningAtDebug)
          printf(
              "[Warning] Cannot hit any more scintillator counters! The velocity in the z-direction is zero or "
              "negative!\n"
          );

        return false;
      }

      // The trajectory is a straight line (x, y, z) = (x0, y0, z0) + t * (vx, vy, vz)
      // The hit time is calculated by the formula: t = (z - z0) / vz
      deltaTime             = (scintillatorCounter.getLocation() - position.Z()) / velocity.Z();
      deltaPropagatedLength = deltaTime * velocity.R();
      hitPosition           = particle->getPosition() + deltaTime * velocity;
      particle->setPosition(hitPosition);

      return true;
    };

    auto spiralPropagate = [&]() -> bool {
      if (1 - this->dB.Unit().Y() > 1e-10) {
        printf(
            "[Error] The magnetic field is not in the y direction! Only supports the magnetic field in the y direction "
            "at the moment!\n"
        );

        exit(EXIT_FAILURE);
      }

      // X-Z plane: (x, z) = (x_c, z_c) + r * (cos(omega * t + theta0), sin(omega * t + theta0))
      // (x_c, z_c) is the center of the circle, r is the radius of the circle, omega is the angular velocity
      // theta0 is associated with the initial position

      double deltaPropagatedLengthXZ = 0, hitPositionX = 0, hitPositionZ = 0, velocityX = 0, velocityZ = 0;
      {
        const ROOT::Math::XYVector velocity0(velocity.X(), velocity.Z()); // initial velocity in the X-Z plane
        if (velocity0.R() < 1e-10) {
          return linearPropagate(); // The velocity is too low, the particle will move in a straight line
        }

        const ROOT::Math::XYPoint position0(position.X(), position.Z()); // initial position in the X-Z plane
        const double              B = this->dB.Y(); // magnetic field (since the magnetic field is in the y direction)
        const double              y = scintillatorCounter.getLocation();

        // F = qv x B = ma, a is always towards the center of the circle
        // v = (vx, vy, 0), B = (0, 0, B) -> qv x B = q(-vy, vx, 0)B
        // so q(vy, -vx, 0) is the direction from the center to the particle
        constexpr double conversionFactor__MeV_c_e_T2cm =
            1e6 / TMath::C() * 1e2; // conversion factor from MeV / c / (e * T) to cm
        constexpr double conversionFactor__cm_ns2c = 1e9 / TMath::Ccgs(); // conversion factor from cm/ns to c
        const double     radius = particle->getMass() * (velocity0.R() * conversionFactor__cm_ns2c) / (B * charge)
                            * conversionFactor__MeV_c_e_T2cm;
        ROOT::Math::XYVector cyclotronDirection(-velocity0.Y(), velocity0.X());
        cyclotronDirection               = cyclotronDirection.Unit();
        const ROOT::Math::XYPoint center = position0 + radius * cyclotronDirection; // radius * cyclotronDirection is the direction from the particle to the center

        const double theta0 = TMath::ATan2(position0.Y() - center.Y(), position0.X() - center.X());
        const double omega  = velocity0.R() / radius;

        double angle = TMath::ASin((y - center.Y()) / radius) - theta0;

        if (std::isnan(angle) || angle < -1e-10) {
          if (Config::enableWarning) {
            printf("[Warning] Cannot hit any more scintillator counters! The velocity might be too low!\n");
          }

          return false;
        }

        deltaTime               = angle / omega;
        deltaPropagatedLengthXZ = radius * angle;
        hitPositionX            = center.X() + radius * TMath::Cos(omega * deltaTime + theta0);
        hitPositionZ            = y;
        velocityX               = -radius * omega * TMath::Sin(omega * deltaTime + theta0);
        velocityZ               = radius * omega * TMath::Cos(omega * deltaTime + theta0);
      }

      // Y direction: y = y0 + vy * t
      const double deltaPropagatedLengthY = velocity.Y() * deltaTime;
      const double hitPositionY           = position.Y() + deltaPropagatedLengthY;

      deltaPropagatedLength = TMath::Sqrt(
          deltaPropagatedLengthXZ * deltaPropagatedLengthXZ + deltaPropagatedLengthY * deltaPropagatedLengthY
      ); // NOTE: this method is ONLY valid for uniform motion

      hitPosition = ROOT::Math::XYZPoint(hitPositionX, hitPositionY, hitPositionZ);
      particle->setPosition(hitPosition);
      particle->setDirection(ROOT::Math::XYZVector(velocityX, velocity.Y(), velocityZ));

      return true;
    };

    if (!charge || this->dB.R() < 1e-10 || position.Z() < initialZ) {
      if (!linearPropagate())
        break; // Magnetic field influence is negligible, the particle will move in a straight line
    } else {
      if (!spiralPropagate()) break; // Magnetic field is not zero, the particle will move in a spiral
    }

    double particleEnergyLoss = 0;
    if (enableEnergyLoss) {
      if (measuresData) particleEnergyLoss = measuresData->at(hitsData->size()).hitEnergyLoss;
      else if (Config::useBetheBloch)
        particleEnergyLoss = scintillatorCounter.energyLoss(*particle);
      else if (Config::enableEnergyLossFluctuation)
        if (enableEnergyLossFluctuation) {
          const double Landau_xi              = scintillatorCounter.LandauMostProbableEnergyLoss_xi(*particle);
          const double mostProbableEnergyLoss = scintillatorCounter.LandauMostProbableEnergyLoss(Landau_xi, *particle);

          if (Config::useLandau) particleEnergyLoss = Config::random->Landau(mostProbableEnergyLoss, 4.018 * Landau_xi);
          else
            particleEnergyLoss = Config::random->Gaus(mostProbableEnergyLoss, 4.018 * Landau_xi);

        } else
          particleEnergyLoss = scintillatorCounter.LandauMostProbableEnergyLoss(*particle);
      else
        particleEnergyLoss = scintillatorCounter.LandauMostProbableEnergyLoss(*particle);

      particle->setEnergy(particle->getEnergy() - particleEnergyLoss);
      if (Config::enableDebug)
        printf(
            "[Info] Energy loss: %f MeV, Energy: %f MeV, Velocity: %f c\n",
            particleEnergyLoss,
            particle->getEnergy(),
            particle->getVelocity().R()
        );
    }

    time += deltaTime, propagationLength += deltaPropagatedLength;
    hitsData->push_back(time, propagationLength, particleEnergyLoss, hitPosition);
  }

  return hitsData;
}

void Detector::plotDeltaTime(const Particle &particle, const std::string &fileName) const {
  const HitsData *hitsDataWithEnergyLoss    = this->particleHitData(particle, true);
  const HitsData *hitsDataWithoutEnergyLoss = this->particleHitData(particle, false);

  this->plotDeltaTime(*hitsDataWithEnergyLoss, *hitsDataWithoutEnergyLoss, fileName);
}

void Detector::plotDeltaTime(
    const HitsData    &hitsDataWithEnergyLoss,
    const HitsData    &hitsDataWithoutEnergyLoss,
    const std::string &fileName
) const {
  TCanvas *Detector_plotDeltaTimeCanvas =
      new TCanvas("Detector_plotDeltaTimeCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
  Detector_plotDeltaTimeCanvas->SetGrid();
  Detector_plotDeltaTimeCanvas->cd();

  TGraph *graph = new TGraph(this->dScintillatorCounters.size());
  graph->SetTitle(";L [cm];#Deltat [ns]");

  for (size_t i = 0; i < this->dScintillatorCounters.size(); ++i) try {
      graph->SetPoint(
          i,
          hitsDataWithEnergyLoss.at(i).hitLength,
          hitsDataWithEnergyLoss.at(i).hitTime - hitsDataWithoutEnergyLoss.at(i).hitTime
      );
    } catch (const std::out_of_range &e) {
      if (Config::enableWarning)
        printf(
            "[Warning] Seems like the particle has not hit the scintillator counter since counter %ld! Only plotting "
            "the hitted scintillator counters!\n",
            i
        );

      break;
    }

  graph->Draw("AL");

  Detector_plotDeltaTimeCanvas->SaveAs(fileName.c_str());

  // Zoom in (Use "{}" to limit the scope of variables)
  {
    double *x = graph->GetX(), *y = graph->GetY();
    int     n  = graph->GetN();
    double  dx = x[n - 1] - x[0], dy = y[n - 1] - y[0];

    std::string baseFileName = fileName.substr(0, fileName.find_last_of('.'));

    if (n > 1) {
      graph->GetXaxis()->SetRangeUser(x[1] - 1e-2 * dx, x[1] + 1e-2 * dx);
      graph->GetYaxis()->SetRangeUser(y[1] - 1e-2 * dy, y[1] + 1e-2 * dy);
      Detector_plotDeltaTimeCanvas->SaveAs(TString::Format("%s_zoom_1.png", baseFileName.c_str()));
    }

    if (n > 2) {
      graph->GetXaxis()->SetRangeUser(x[2] - 1e-2 * dx, x[2] + 1e-2 * dx);
      graph->GetYaxis()->SetRangeUser(y[2] - 1e-2 * dy, y[2] + 1e-2 * dy);
      Detector_plotDeltaTimeCanvas->SaveAs(TString::Format("%s_zoom_2.png", baseFileName.c_str()));
    }
  }

  delete graph;
  delete Detector_plotDeltaTimeCanvas;
}

HitsData *Detector::measure(const HitsData &hitsData) const {
  HitsData *measuresData = MemoryPool::getDetector_detect_detectedsData(this->dScintillatorCounters.size());
  *measuresData          = hitsData;

  for (size_t i = 0; i < measuresData->size(); ++i) {
    if (Config::enableTimeResolution)
      measuresData->setHitTime(
          i,
          Config::random->Gaus(measuresData->at(i).hitTime, this->dScintillatorCounters.at(i).getTimeResolution())
      );

    if (!Config::useRealEnergyLoss2Rec)
      measuresData->setHitEnergyLoss(i, Config::random->Gaus(measuresData->at(i).hitEnergyLoss, 20));
  }

  return measuresData;
}

double Detector::reconstructUsingLinearMethod(const HitsData &measuresData) const {
  const int n = measuresData.size();

  // Fit t = kL + b using least squares method
  // we only need the slope k
  // k = [n*sum(L*t) - sum(L)*sum(t)] / [n*sum(L^2) - (sum(L))^2]
  double sumL  = 0.0; // sum(L)
  double sumT  = 0.0; // sum(t)
  double sumLT = 0.0; // sum(L*t)
  double sumL2 = 0.0; // sum(L^2)

  for (int i = 0; i < n; ++i) {
    const double L  = measuresData.at(i).hitLength;
    const double t  = measuresData.at(i).hitTime;
    sumL           += L;
    sumT           += t;
    sumLT          += L * t;
    sumL2          += L * L;
  }

  const double denominator = n * sumL2 - sumL * sumL;
  if (std::abs(denominator) < 1e-10) {
    if (Config::enableWarning) printf("[Warning] The particle is too slow!");
    return 0;
  }

  const double k = (n * sumLT - sumL * sumT) / denominator;

  constexpr double conversionFactor__ns_cm2_c = TMath::Ccgs() * 1e-9; // conversion factor from ns/cm to 1/c
  const double     betaReciprocal             = k * conversionFactor__ns_cm2_c;

  if (Config::enableDebug) printf("[Info] The reconstructed 1/beta using the linear method: %f\n", betaReciprocal);

  return betaReciprocal;
}

void Detector::plotReconstructDataUsingLinearMethod(const Particle &particle, const std::string &fileName) const {
  const HitsData *hitsData     = this->particleHitData(particle, true, false);
  const HitsData *measuresData = this->measure(*hitsData);
  const int       n            = measuresData->size();

  const std::vector<double> &hitsTime        = hitsData->getHitsTime();
  const std::vector<double> &hitsLengths     = hitsData->getHitsLength();
  const std::vector<double> &measuresTime    = measuresData->getHitsTime();
  const std::vector<double> &measuresLengths = measuresData->getHitsLength();

  TCanvas *Detector_plotReconstructDataCanvas =
      new TCanvas("Detector_plotReconstructDataCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
  Detector_plotReconstructDataCanvas->SetGrid();
  Detector_plotReconstructDataCanvas->cd();

  TLegend *legend = new TLegend(0.2, 0.75, 0.5, 0.85);
  legend->SetBorderSize(kNone);

  TGraph *graphRealData = new TGraph(n, &hitsLengths[0], &hitsTime[0]);
  graphRealData->SetMarkerStyle(kFullTriangleUp);
  graphRealData->SetMarkerSize(3);
  graphRealData->SetTitle(";L [cm];t [ns]");
  legend->AddEntry(graphRealData, "Real", "p");
  graphRealData->Draw("AP");

  TGraph *graphReconstructData = new TGraph(n, &measuresLengths[0], &measuresTime[0]);
  graphReconstructData->SetMarkerStyle(kFullCircle);
  graphReconstructData->SetMarkerColor(kRed);
  graphReconstructData->SetMarkerSize(2);
  legend->AddEntry(graphReconstructData, "Detected", "p");
  graphReconstructData->Draw("PSAME");

  TF1 *f1 = new TF1("f1", "[0] * x", 0, measuresLengths.back());
  graphReconstructData->Fit(f1, "Q");
  const int    exponent = TMath::FloorNint(TMath::Log10(f1->GetParameter(0)));
  const double mantissa = f1->GetParameter(0) / TMath::Power(10, exponent);
  legend->AddEntry(f1, Form("y = %.3f#times10^{%i} x (fit line)", mantissa, exponent));

  legend->Draw();

  if (Config::enableDebug) {
    constexpr double conversionFactor__ns_cm2_c = TMath::Ccgs() * 1e-9; // conversion factor from ns/cm to 1/c
    const double     betaReciprocal             = f1->GetParameter(0) * conversionFactor__ns_cm2_c;
    printf("[Info] The reconstructed 1 / beta using the linear method: %f\n", betaReciprocal);
    printf("[Info] The real 1 / beta: %f\n", 1 / particle.getBeta());
  }

  Detector_plotReconstructDataCanvas->SaveAs(fileName.c_str());

  delete graphRealData;
  delete graphReconstructData;
  delete f1;
  delete legend;
  delete Detector_plotReconstructDataCanvas;
}

std::pair<double, double> Detector::distributionOfReconstruction(
    const Particle    &particle,
    const int          nReconstructions,
    const bool         enableLinearMethod,
    const bool         enablePlot,
    const std::string &fileName
) const {
  const double betaReciprocalReal = 1 / particle.getBeta();

  std::vector<double> results(nReconstructions, 0);

  for (int i = 0; i < nReconstructions; ++i) {
    HitsData *hitsData =
        MemoryPool::getDetector_distributionOfReconstruction_hitsData(this->dScintillatorCounters.size());
    *hitsData = *this->particleHitData(particle, true);

    if (hitsData->size() != this->dScintillatorCounters.size()) {
      if (Config::enableWarning)
        printf(
            "[Warning] Cancel reconstruction No.%d since the particle has not hit all the scintillator counters!\n",
            i
        );

      results[i] = std::nan("");
      continue;
    }

    double betaReciprocalReconstruction;
    if (enableLinearMethod)
      betaReciprocalReconstruction = this->reconstructUsingLinearMethod(*this->measure(*hitsData));
    else
      betaReciprocalReconstruction = this->reconstructUsingNonLinearMethod(particle, *this->measure(*hitsData));

    results[i] = betaReciprocalReal - betaReciprocalReconstruction;
  }

  double deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation;
  {
    const HitsData *hitsData = this->particleHitData(particle, true, false);
    size_t          n        = hitsData->size();

    // adjust the beta of the particle to make sure it can hit all the scintillator counters
    if (n != this->dScintillatorCounters.size()) {
      if (Config::enableWarning)
        printf(
            "[Warning] Cannot reconstruct 1 / beta with zero resolution and zero energy loss fluctuation since the "
            "particle has not hit all the scintillator counters at beta = %f!\n",
            particle.getBeta()
        );

      printf("[Info] Try to find the minimum beta of the particle that can hit all the scintillator counters ...\n");
      double   betaMin = particle.getBeta(), betaMax = 0.99;
      Particle particleTemp(particle);
      while (abs(betaMax - betaMin) > 1e-3) {
        double beta = (betaMin + betaMax) / 2;
        particleTemp.setBeta(beta);
        hitsData = this->particleHitData(particleTemp, true, false);
        n        = hitsData->size();
        if (n == this->dScintillatorCounters.size()) betaMax = beta;
        else
          betaMin = beta;
      }
      if (betaMax >= 0.988) {
        printf("[Error] The minimum beta is too close to 1! The particle cannot hit all the scintillator counters!\n");
        exit(EXIT_FAILURE);
      } else
        printf("[Info] The minimum beta of the particle that can hit all the scintillator counters: %f\n", betaMax);

      if (Config::enableWarning)
        printf(
            "[Warning] Use %f as the beta of the particle to reconstruct the 1 / beta with zero resolution and zero "
            "energy loss fluctuation!\n",
            betaMax
        );

      particleTemp.setBeta(betaMax);
      hitsData = this->particleHitData(particleTemp, true, false);
      n        = hitsData->size();
    }

    if (enableLinearMethod)
      deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation =
          betaReciprocalReal - this->reconstructUsingLinearMethod(*hitsData);
    else
      deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation =
          betaReciprocalReal - this->reconstructUsingNonLinearMethod(particle, *hitsData);
  }

  // Sometimes, the histogram seems not to be deleted correctly, so we use a counter to avoid the warning of ROOT
  static int histCounter = 0;

  TH1F *histogram = new TH1F(
      Form("#Delta(1/#beta)_%d", histCounter++),
      ";#Delta(1/#beta);",
      1000,
      deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation - 1,
      deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation + 1
  );
  for (const auto &result : results)
    if (!std::isnan(result)) histogram->Fill(result);

  TCanvas *Detector_plotDistributionOfReconstructionCanvas;
  if (enablePlot) {
    Detector_plotDistributionOfReconstructionCanvas =
        new TCanvas("Detector_plotDistributionOfReconstructionCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotDistributionOfReconstructionCanvas->cd();

    TFitResultPtr fitResult = histogram->Fit("gaus", "QS");
    histogram->GetFunction("gaus")->SetNpx(1000);

    double mean = fitResult->Value(1), sigma = fitResult->Value(2);

    histogram->GetXaxis()->SetRangeUser(
        std::max(mean - 5 * sigma, histogram->GetXaxis()->GetXmin()),
        std::min(mean + 5 * sigma, histogram->GetXaxis()->GetXmax())
    );

    histogram->Draw();

    TLine *line = new TLine(
        deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation,
        0,
        deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation,
        histogram->GetMaximum()
    );
    line->Draw("same");

    TLegend *legend = new TLegend(0.62, 0.67, 0.87, 0.87);
    legend->SetBorderSize(kNone);
    legend->AddEntry("", Form("#chi^{2} / NDF = %.4g / %u", fitResult->Chi2(), fitResult->Ndf()), "");
    legend->AddEntry("", Form("#mu = %.4g", mean), "");
    legend->AddEntry("", Form("#sigma = %.4g", sigma), "");
    legend->AddEntry("", Form("Entries: %d", (int)histogram->GetEntries()), "");
    legend->AddEntry(
        line,
        Form("#Delta(1/#beta)_{real} = %.4g", deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation),
        "l"
    );
    legend->Draw();

    if (Config::enableDebug) {
      printf("[Info] The mean of the distribution of the difference between real and reconstructed 1/beta: %f\n", mean);
      printf(
          "[Info] The difference between real and reconstructed 1/beta with zero resolution and zero energy loss "
          "fluctuation: %f\n",
          deltaBetaReciprocalWithZeroResolutionAndZeroEnergyLossFluctuation
      );
    }

    Detector_plotDistributionOfReconstructionCanvas->SaveAs(fileName.c_str());

    delete line;
    delete legend;
    delete Detector_plotDistributionOfReconstructionCanvas;
  }

  std::pair<double, double> meanAndStandardDeviation = {histogram->GetMean(), histogram->GetRMS()};

  delete histogram;

  return meanAndStandardDeviation;
}

TGraphErrors *Detector::deltaBetaReciprocal(
    Particle           particle,
    const double       betaMin,
    const double       betaMax,
    const int          nPoints,
    const bool         enableLinearMethod,
    const bool         enablePlot,
    const std::string &fileName
) const {
  TCanvas *Detector_plotDeltaBetaReciprocalCanvas = nullptr;
  if (enablePlot) {
    Detector_plotDeltaBetaReciprocalCanvas =
        new TCanvas("Detector_plotDeltaBetaReciprocalCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
    Detector_plotDeltaBetaReciprocalCanvas->SetGrid();
    Detector_plotDeltaBetaReciprocalCanvas->cd();
  }

  TGraphErrors *graphErrors = new TGraphErrors(nPoints + 1);
  graphErrors->SetTitle(";#beta;#Delta(1/#beta)");

  const double step = (betaMax - betaMin) / nPoints;

  std::vector<std::pair<double, double>> deltaBetaReciprocals;
  deltaBetaReciprocals.reserve(nPoints + 1);

  if (Config::enableMultiThreading) {
    std::vector<std::future<std::pair<double, double>>> futures;
    futures.reserve(nPoints + 1);
    for (int i = 0; i <= nPoints; ++i) {
      const double beta = betaMin + i * step;
      particle.setBeta(beta);

      futures.push_back(
          ThreadPool::getInstance().enqueue([this, particle, enableLinearMethod]() -> std::pair<double, double> {
            return this->distributionOfReconstruction(particle, 10000, enableLinearMethod, false);
          })
      );
    }

    for (auto &future : futures) deltaBetaReciprocals.push_back(future.get());
  } else
    for (int i = 0; i <= nPoints; ++i) {
      const double beta = betaMin + i * step;
      particle.setBeta(beta);

      deltaBetaReciprocals.push_back(this->distributionOfReconstruction(particle, 10000, enableLinearMethod, false));
    }

  for (int i = 0; i <= nPoints; ++i) {
    std::pair<double, double> &deltaBetaReciprocal = deltaBetaReciprocals.at(i);
    const double               beta                = betaMin + i * step;
    graphErrors->SetPoint(i, beta, deltaBetaReciprocal.first);
    graphErrors->SetPointError(i, 0, deltaBetaReciprocal.second);
  }

  if (enablePlot) {
    graphErrors->Draw("AL");
    Detector_plotDeltaBetaReciprocalCanvas->SaveAs(fileName.c_str());

    delete Detector_plotDeltaBetaReciprocalCanvas;
  }

  return graphErrors;
}

double
Detector::reconstructUsingNonLinearMethod(const Particle &particleOrignal, const HitsData &measuresDataOrignal) const {
  Particle *particle = MemoryPool::getDetector_reconstructUsingNonLinearMethod_particle();
  *particle          = particleOrignal;

  HitsData *measuresData =
      MemoryPool::getDetector_reconstructUsingNonLinearMethod_measuresData(this->getScintillatorCounters().size());
  *measuresData = measuresDataOrignal;

  double initialBetaReciprocal = this->reconstructUsingLinearMethod(*measuresData);

  auto chi2Function = [&](const double *parameters) {
    particle->setBeta(1 / parameters[0]);

    HitsData *reconstructedHitsData;
    if (Config::useEnergyLoss2Rec) reconstructedHitsData = this->particleHitData(*particle, true, false, measuresData);
    else
      reconstructedHitsData = this->particleHitData(*particle, true, false);

    if (reconstructedHitsData->size() != this->getScintillatorCounters().size()) {
      if (Config::enableWarningAll)
        printf(
            "[Warning] [@Reconstruction] The particle has not hit all the scintillator counters as beta = %f!\n",
            1 / parameters[0]
        );

      return 1e10 * parameters[0];
    }

    double chi2 = 0;
    for (size_t i = 0; i < reconstructedHitsData->size(); ++i) {
      const double timeDiff        = reconstructedHitsData->at(i).hitTime - measuresData->at(i).hitTime;
      const double timeResolution  = this->getScintillatorCounters().at(i).getTimeResolution();
      chi2                        += TMath::Power(timeDiff, 2) / TMath::Power(timeResolution, 2);
    }
    return chi2;
  };

  ROOT::Math::Minimizer *minimizer = MemoryPool::getDetector_reconstructUsingNonLinearMethod_minimizer();

  ROOT::Math::Functor functor(chi2Function, 1);
  minimizer->SetFunction(functor);

  initialBetaReciprocal =
      std::clamp(initialBetaReciprocal, 1 + 1e-4, 10 - 1e-4); // restrict the range of betaReciprocal between 1 and 10
  minimizer->SetLimitedVariable(0, "betaReciprocal", initialBetaReciprocal, 1e-5, 1, 10);

  if (Config::enableDebug) minimizer->SetPrintLevel(1);

  minimizer->Minimize();
  const double betaReciprocal = minimizer->X()[0];

  if (Config::enableDebug) {
    printf("[Info] The reconstructed 1/beta using the linear method: %f\n", initialBetaReciprocal);
    printf("[Info] The reconstructed 1/beta using the non-linear method: %f\n", betaReciprocal);
  }

  return betaReciprocal;
}

void Detector::plotParticleTrajectory(Particle particle, const std::string &fileName) const {
  TCanvas *canvas = new TCanvas("Detector_plotParticleTrajectoryCanvas", "", 3508, 2480); // A4 size in pixels(300 dpi)
  canvas->Divide(2, 2);

  const HitsData *hitsData = this->particleHitData(particle, true, false);
  const int       nPoints  = hitsData->size();
  if (nPoints < 2) {
    if (Config::enableWarning) { printf("[Warning] Not enough points to plot trajectory (nPoints = %d)\n", nPoints); }
    return;
  }

  TMultiGraph *multiGraphXY = new TMultiGraph();
  TMultiGraph *multiGraphYZ = new TMultiGraph();
  TMultiGraph *multiGraphXZ = new TMultiGraph();

  const int totalSamples = (nPoints - 1) * 100;
  TGraph2D *graph3D      = new TGraph2D(totalSamples);
  int       pointIndex   = 0;

  for (int i = 0; i < nPoints - 1; ++i) {
    const ROOT::Math::XYZPoint &pos1 = hitsData->at(i).hitPosition;
    const ROOT::Math::XYZPoint &pos2 = hitsData->at(i + 1).hitPosition;

    const double energyLoss = hitsData->at(i).hitEnergyLoss;
    particle.setEnergy(particle.getEnergy() - energyLoss);
    if (Config::enableDebug) {
      printf("[Info] Energy loss at point %d: %f MeV, Remaining energy: %f MeV\n", i, energyLoss, particle.getEnergy());
    }

    const int nSamples = 100;
    double   *xTrack   = new double[nSamples];
    double   *yTrack   = new double[nSamples];
    double   *zTrack   = new double[nSamples];

    constexpr double            conversionFactor__c2cm_ns = TMath::Ccgs() * 1e-9; // conversion factor from c to cm/ns
    const ROOT::Math::XYZVector velocity                  = particle.getVelocity() * conversionFactor__c2cm_ns;
    const int                   charge                    = particle.getCharge();
    const bool                  hasMagneticField          = this->dB.R() >= 1e-10;

    if (!charge || !hasMagneticField) {
      for (int j = 0; j < nSamples; ++j) {
        const double t = j / (nSamples - 1.0);
        xTrack[j]      = pos1.X() + (pos2.X() - pos1.X()) * t;
        yTrack[j]      = pos1.Y() + (pos2.Y() - pos1.Y()) * t;
        zTrack[j]      = pos1.Z() + (pos2.Z() - pos1.Z()) * t;
      }
    } else {
      const double B = this->dB.Y();

      const ROOT::Math::XYVector velocity0(velocity.X(), velocity.Z());
      if (velocity0.R() < 1e-10) {
        for (int j = 0; j < nSamples; ++j) {
          const double t = j / (nSamples - 1.0);
          xTrack[j]      = pos1.X() + (pos2.X() - pos1.X()) * t;
          yTrack[j]      = pos1.Y() + (pos2.Y() - pos1.Y()) * t;
          zTrack[j]      = pos1.Z() + (pos2.Z() - pos1.Z()) * t;
        }
      } else {
        constexpr double conversionFactor__MeV_c_e_T2cm = 1e6 / TMath::C() * 1e2; // conversion factor from MeV/c to cm
        constexpr double conversionFactor__cm_ns2c      = 1e9 / TMath::Ccgs();    // conversion factor from cm/ns to c
        const double     radius = particle.getMass() * (velocity0.R() * conversionFactor__cm_ns2c) / (B * abs(charge))
                            * conversionFactor__MeV_c_e_T2cm;

        ROOT::Math::XYVector cyclotronDirection(velocity0.Y(), -velocity0.X());
        cyclotronDirection = cyclotronDirection.Unit();
        const ROOT::Math::XYPoint center(
            pos1.X() - radius * cyclotronDirection.X(),
            pos1.Z() - radius * cyclotronDirection.Y()
        );

        const double theta0    = TMath::ATan2(pos1.Z() - center.Y(), pos1.X() - center.X());
        const double omega     = velocity0.R() / radius;
        const double deltaTime = (TMath::ASin((pos2.Z() - center.Y()) / radius) - theta0) / omega;

        for (int j = 0; j < nSamples; ++j) {
          const double t = j * deltaTime / (nSamples - 1.0);
          xTrack[j]      = center.X() + radius * TMath::Cos(omega * t + theta0);
          yTrack[j]      = pos1.Y() + velocity.Y() * t;
          zTrack[j]      = center.Y() + radius * TMath::Sin(omega * t + theta0);
        }

        const double velocityX = -radius * omega * TMath::Sin(omega * deltaTime + theta0);
        const double velocityZ = radius * omega * TMath::Cos(omega * deltaTime + theta0);
        particle.setVelocity(ROOT::Math::XYZVector(velocityX, velocity.Y(), velocityZ) / conversionFactor__c2cm_ns);
      }
    }

    TGraph *segmentXY = new TGraph(nSamples, xTrack, yTrack);
    TGraph *segmentYZ = new TGraph(nSamples, zTrack, yTrack);
    TGraph *segmentXZ = new TGraph(nSamples, zTrack, xTrack);

    segmentXY->SetLineColor(kBlue);
    segmentYZ->SetLineColor(kBlue);
    segmentXZ->SetLineColor(kBlue);

    multiGraphXY->Add(segmentXY);
    multiGraphYZ->Add(segmentYZ);
    multiGraphXZ->Add(segmentXZ);

    for (int j = 0; j < nSamples; ++j) { graph3D->SetPoint(pointIndex++, xTrack[j], yTrack[j], zTrack[j]); }

    delete[] xTrack;
    delete[] yTrack;
    delete[] zTrack;
  }

  canvas->cd(1);
  multiGraphXY->SetTitle("XY Plane View (Top);X [cm];Y [cm]");
  multiGraphXY->Draw("AL");
  gPad->SetGrid();

  canvas->cd(2);
  multiGraphYZ->SetTitle("YZ Plane View (Side);Z [cm];Y [cm]");
  multiGraphYZ->Draw("AL");
  gPad->SetGrid();

  canvas->cd(3);
  multiGraphXZ->SetTitle("XZ Plane View (Front);Z [cm];X [cm]");
  multiGraphXZ->Draw("AL");
  gPad->SetGrid();

  canvas->cd(4);
  graph3D->SetTitle("3D View;X [cm];Y [cm];Z [cm]");
  graph3D->SetLineColor(kBlue);
  graph3D->SetLineWidth(2);
  graph3D->Draw("LINE");

  gPad->SetTheta(30);
  gPad->SetPhi(30);
  gPad->SetGrid();

  canvas->SaveAs(fileName.c_str());

  delete multiGraphXY;
  delete multiGraphYZ;
  delete multiGraphXZ;
  delete graph3D;
  delete canvas;
}
