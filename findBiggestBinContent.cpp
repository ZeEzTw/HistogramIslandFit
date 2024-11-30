{
    const int ndet = 32;
    const int nE = 2;
    const double Eg[nE] = {1173, 1332};

    // Declare arrays to store fit parameters
    double fpar0[ndet][ndet][nE];
    double fpar1[ndet][ndet][nE];
    double fparerr0[ndet][ndet][nE];
    double fparerr1[ndet][ndet][nE];

    TCanvas *canvas = new TCanvas("c", "Histograms", 800, 600);
    TFile *outFile = new TFile("output/histograms_with_lines_and_slopes.root", "RECREATE");

    std::map<std::pair<int, int>, double> errorMap;
    std::ifstream infile("data/CL32_matrix/CL32_matrix.dat");
    int i, j;
    double errorValue;

    while (infile >> i >> j >> errorValue)
    {
        errorMap[std::make_pair(i, j)] = errorValue;
    }

    std::ofstream resultFile("output/results.txt");

    for (int i = 0; i < ndet; i++)
    {
        resultFile <<1.0 << "\n";

        for (int j = i + 1; j < ndet; j++)
        {
            TFile *fEE = new TFile(Form("data/hEE_%i_%i.root", i, j));
            if (!fEE || fEE->IsZombie())
            {
                delete fEE;
                continue;
            }
            TH2D *hEE;
            fEE->GetObject(Form("hEE%i %i", i, j), hEE);
            if (!hEE)
            {
                delete fEE;
                continue;
            }

            hEE->GetXaxis()->SetRangeUser(0, 1500);
            hEE->GetYaxis()->SetRangeUser(0, 1500);
            hEE->Draw("COLZ");

            auto key1 = std::make_pair(i + 1, j + 1);
            auto key2 = std::make_pair(j + 1, i + 1);

            if (errorMap.find(key1) != errorMap.end() && errorMap.find(key2) != errorMap.end())
            {
                double errorValue1 = errorMap[key1];
                double xStart1 = 1173 + errorValue1 * 1332 - 100;
                double yStart1 = 1332 + errorValue1 * 1173 - 200;
                double xCoord1 = xStart1, yCoord1 = yStart1;
                cout << xStart1 << " " << yStart1 << endl;
                double largestSum1 = 0;
                for (int offsetX = 0; offsetX < 150; ++offsetX)
                {
                    for (int offsetY = 0; offsetY < 250; ++offsetY)
                    {
                        double currentSum = 0;
                        for (int binX = xStart1 + offsetX; binX < xStart1 + offsetX + 7; ++binX)
                        {
                            for (int binY = yStart1 + offsetY; binY < yStart1 + offsetY + 7; ++binY)
                            {
                                currentSum += hEE->GetBinContent(hEE->GetXaxis()->FindBin(binX), hEE->GetYaxis()->FindBin(binY));
                            }
                        }
                        if (currentSum > largestSum1)
                        {
                            largestSum1 = currentSum;
                            xCoord1 = xStart1 + offsetX;
                            yCoord1 = yStart1 + offsetY;
                        }
                    }
                }
                
                TF2 *fitFunc1 = new TF2("fitFunc1", "[0]*exp(-0.5*((x-[1])/[2])^2 - 0.5*((y-[3])/[4])^2)", 
                                        xCoord1 - 15, xCoord1 + 15, yCoord1 - 15, yCoord1 + 15);
                
                double maxContent1 = hEE->GetBinContent(hEE->GetXaxis()->FindBin(xCoord1), 
                                                      hEE->GetYaxis()->FindBin(yCoord1));
                
                fitFunc1->SetParameters(maxContent1, xCoord1, 2.5, yCoord1, 2.5);
                // Set parameter bounds
                fitFunc1->SetParLimits(0, 0.5*maxContent1, 2*maxContent1); // amplitude
                fitFunc1->SetParLimits(1, xCoord1-5, xCoord1+5);  // mean x
                fitFunc1->SetParLimits(2, 1.0, 5.0);              // sigma x
                fitFunc1->SetParLimits(3, yCoord1-5, yCoord1+5);  // mean y
                fitFunc1->SetParLimits(4, 1.0, 5.0);              // sigma y
                
                hEE->Fit(fitFunc1, "QR");
                
                double fitX1 = fitFunc1->GetParameter(1);
                double fitY1 = fitFunc1->GetParameter(3);
                cout<<"Fit X1: "<<fitX1<<" Fit Y1: "<<fitY1<<endl;

                xCoord1 = fitX1;
                yCoord1 = fitY1;

                double errorValue2 = errorMap[key2];
                double xStart2 = 1332 + errorValue2 * 1173 - 100;
                double yStart2 = 1173 + errorValue2 * 1332 - 200;
                double xCoord2 = xStart2, yCoord2 = yStart2;

                double largestSum2 = 0;
                for (int offsetX = 0; offsetX < 150; ++offsetX)
                {
                    for (int offsetY = 0; offsetY < 250; ++offsetY)
                    {
                        double currentSum = 0;
                        for (int binX = xStart2 + offsetX; binX < xStart2 + offsetX + 7; ++binX)
                        {
                            for (int binY = yStart2 + offsetY; binY < yStart2 + offsetY + 7; ++binY)
                            {
                                currentSum += hEE->GetBinContent(hEE->GetXaxis()->FindBin(binX), hEE->GetYaxis()->FindBin(binY));
                            }
                        }
                        if (currentSum > largestSum2)
                        {
                            largestSum2 = currentSum;
                            xCoord2 = xStart2 + offsetX;
                            yCoord2 = yStart2 + offsetY;
                        }
                    }
                }

                TF2 *fitFunc2 = new TF2("fitFunc2", "[0]*exp(-0.5*((x-[1])/[2])^2 - 0.5*((y-[3])/[4])^2)", 
                                        xCoord2 - 15, xCoord2 + 15, yCoord2 - 15, yCoord2 + 15);
                
                double maxContent2 = hEE->GetBinContent(hEE->GetXaxis()->FindBin(xCoord2), 
                                                      hEE->GetYaxis()->FindBin(yCoord2));
                
                fitFunc2->SetParameters(maxContent2, xCoord2, 2.5, yCoord2, 2.5);

                fitFunc2->SetParLimits(0, 0.5*maxContent2, 2*maxContent2); // amplitude
                fitFunc2->SetParLimits(1, xCoord2-5, xCoord2+5);  // mean x
                fitFunc2->SetParLimits(2, 1.0, 5.0);              // sigma x
                fitFunc2->SetParLimits(3, yCoord2-5, yCoord2+5);  // mean y
                fitFunc2->SetParLimits(4, 1.0, 5.0);              // sigma y
                
                hEE->Fit(fitFunc2, "QR");
                
                double fitX2 = fitFunc2->GetParameter(1);
                double fitY2 = fitFunc2->GetParameter(3);
                
                xCoord2 = fitX2;
                yCoord2 = fitY2;

                TEllipse *circle1 = new TEllipse(xCoord1, yCoord1, 7, 7);
                circle1->SetLineColor(kRed);
                circle1->SetFillStyle(0);
                circle1->SetLineWidth(2);
                circle1->Draw("same");

                TEllipse *circle2 = new TEllipse(xCoord2, yCoord2, 7, 7);
                circle2->SetLineColor(kGreen);
                circle2->SetFillStyle(0);
                circle2->SetLineWidth(2);
                circle2->Draw("same");

                double xCoords[2] = {xCoord1, xCoord2};
                double yCoords[2] = {yCoord1, yCoord2};
                TLine *line = new TLine(xCoord1, yCoord1, xCoord2, yCoord2);
                line->SetLineColor(kBlue);
                line->SetLineWidth(2);
                line->Draw("same");

                // Adjusted line fitting
                TGraph *adjustmentLine = new TGraph(2, xCoords, yCoords);

                TF1 *f1 = new TF1("f1", "[0]*x + [1]", 0, 1500);
                adjustmentLine->Fit("f1", "Q");

                //adjustmentLine->Draw("same");
                double slope = f1->GetParameter(0);         // Slope
                double intercept = f1->GetParameter(1);     // Intercept
                double slopeError = f1->GetParError(0);     // Slope error
                double interceptError = f1->GetParError(1); // Intercept error

                // Store fit parameters
                for (int k = 0; k < nE; ++k)
                {
                    fpar0[i][j][k] = intercept;
                    fpar1[i][j][k] = slope;
                    fparerr0[i][j][k] = interceptError;
                    fparerr1[i][j][k] = S;
                }

                double delta_ij = 0;
                for (int k = 0; k < nE; k++) {
                    delta_ij += (fpar0[i][j][k]/Eg[k] + fpar1[i][j][k]); // First formula
                }
                delta_ij = delta_ij/nE - 1.0; // Average and subtract 1.0 at the end

                resultFile<<delta_ij << "\n";
                
                // Calculate and write the reverse pair (j,i)
                double delta_ji = 0;
                for (int k = 0; k < nE; k++) {
                    delta_ji += (1.0 - fpar0[i][j][k]/Eg[k])/fpar1[i][j][k];
                }
                delta_ji = delta_ji/nE - 1.0;

                resultFile<< delta_ji << "\n";

                TLatex *latex = new TLatex();
                latex->SetTextSize(0.025);
                latex->DrawLatex((xCoord1 + xCoord2) / 2, (yCoord1 + yCoord2) / 2, Form("Slope = %.4f", slope));

                canvas->Update();
                canvas->SaveAs(Form("output/histogram_with_line_%i_%i.root", i, j));
                outFile->WriteTObject(hEE, Form("histogram_with_line_%i_%i", i, j));

                // resultFile << i + 1 << " " << j + 1 << " " << slope << "\n";
            }

            delete hEE;
            delete fEE;
        }
    }

    resultFile.close();
    outFile->Close();
    return 0;
}
