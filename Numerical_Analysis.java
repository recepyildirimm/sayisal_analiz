class Numerical_Analysis {

    DecimalFormat frm = new DecimalFormat("#.####");

    void ileriSonluFarkHesapla(ArrayList<Double> arrayList, double h) {
        double deltaY;
        ArrayList<Double> deltaYList = new ArrayList<>();
        for (int i = 0; i < arrayList.size(); i++) {
            System.out.println((i + 1) + ".MERTEBEDEN DELTA Y LER:");
            if (deltaYList.size() > 0) {
                for (int j = 0; j < deltaYList.size() - 1; j++) {
                    deltaY = deltaYList.get(j + 1) - deltaYList.get(j);
                    deltaYList.set(j, deltaY);
                    System.out.println("Delta Y:" + (j) + "=" + frm.format(deltaY) + "---");
                }
                deltaYList.remove(deltaYList.size() - 1);
            } else {
                for (int j = 0; j < arrayList.size() - 1; j++) {
                    deltaY = arrayList.get(j + 1) - arrayList.get(j);
                    System.out.println("Delta Y:" + (j) + "=" + frm.format(deltaY) + "---");
                    deltaYList.add(deltaY);

                }

            }
        }
    }

    void geriSonluFarkHesapla(ArrayList<Double> arrayList, double h) {
        double deltaY;
        ArrayList<Double> deltaYList = new ArrayList<>();
        for (int i = 0; i < arrayList.size(); i++) {
            System.out.println((i + 1) + ".MERTEBEDEN DELTA Y LER:");
            if (deltaYList.size() > 0) {
                for (int j = 1; j < deltaYList.size(); j++) {
                    deltaY = (deltaYList.get(j - 1) - deltaYList.get(j));
                    deltaYList.set(j - 1, deltaY);
                    System.out.println("Delta Y:" + (arrayList.size() - j - 1) + "=" + frm.format(deltaY) + "---");
                }
                deltaYList.remove(deltaYList.size() - 1);
            } else {
                for (int j = arrayList.size() - 1; j > 0; j--) {
                    deltaY = (arrayList.get(j) - arrayList.get(j - 1));
                    System.out.println("Delta Y:" + (j) + "=" + frm.format(deltaY) + "---");
                    deltaYList.add(deltaY);

                }

            }
        }
    }

    void merkeziSonluFarkHesapla(ArrayList<Double> arrayList, double h) {
        double deltaY;
        ArrayList<Double> deltaYList = new ArrayList<>();
        for (int i = 0; i < arrayList.size(); i++) {
            System.out.println((i + 1) + ".MERTEBEDEN DELTA Y LER:");
            if (deltaYList.size() > 0) {
                for (int j = 0; j < deltaYList.size() - 1; j++) {
                    deltaY = deltaYList.get(j + 1) - deltaYList.get(j);
                    deltaYList.set(j, deltaY);
                    System.out.println("Delta Y:" + (j + 0.5 + (i * 0.5)) + "=" + frm.format(deltaY) + "---");
                }
                deltaYList.remove(deltaYList.size() - 1);
            } else {
                for (int j = 0; j < arrayList.size() - 1; j++) {
                    deltaY = arrayList.get(j + 1) - arrayList.get(j);
                    System.out.println("Delta Y:" + (j + 0.5 + (i * 0.5)) + "=" + frm.format(deltaY) + "---");
                    deltaYList.add(deltaY);

                }

            }
        }
    }

    void interpolasyonHesapla(double x0, double y0, double deltaY0, double x, double gercekDeger, double h) {
        double s = (x - x0) / h;
        double sonuc = y0 + deltaY0 * s;
        System.out.println("f(" + x + ") = " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void quadratikInterpolasyonHesapla(double x0, double y0, double deltaY0, double deltaY02, double x, double gercekDeger, double h) {
        double s = (x - x0) / h;
        double sonuc = y0 + deltaY0 * s + ((s * (s - 1)) / 2) * deltaY02;
        System.out.println("f(" + x + ") = " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

      void hsizQuadratikInterpolasyonHesapla(double x0, double y0, double x1, double y1, double x2, double y2, double soruX, double gercekDeger) {
        double b0 = y0;
        double b1 = (y1 - y0) / (x1 - x0);
        double b2 = (((y2 - y1) / (x2 - x1)) - ((y1 - y0) / (x1 - x0))) / (x2 - x0);
        double sonuc = b0 + (b1 * (soruX - x0)) + (b2 * (soruX - x0) * (soruX - x1));
        System.out.println("b0 = " + frm.format(b0) + " b1 = " + frm.format(b1) + " b2 = " + frm.format(b2));
        System.out.println("f(" + soruX + ") = " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }
    void newtonInterpolasyonHesapla(ArrayList<Double> xDegerleri, ArrayList<Double> yDegerleri, double soruX, double gercekDeger) {
        ArrayList<Double> deltayDegerleri = new ArrayList<>();
        ArrayList<Double> bDegerleri = new ArrayList<>();
        for (int i = 0; i < xDegerleri.size(); i++) {
            if (deltayDegerleri.size() > 0) {
                for (int j = 0; j < deltayDegerleri.size() - 1; j++) {
                    double nYDegeri = (deltayDegerleri.get(j + 1) - deltayDegerleri.get(j)) / (xDegerleri.get(j + i) - xDegerleri.get(j));
                    deltayDegerleri.set(j, nYDegeri);
                    System.out.println((i + 1) + "indexli ().+" + (j + 1) + "f fonksiyonu = " + frm.format(nYDegeri));
                }
                System.out.println("b" + i + " = " + frm.format(deltayDegerleri.get(0)));
                bDegerleri.add(deltayDegerleri.get(0));
                deltayDegerleri.remove(deltayDegerleri.size() - 1);
            } else {
                for (int j = 0; j < yDegerleri.size(); j++) {
                    deltayDegerleri.add(yDegerleri.get(j));
                }
                System.out.println("b0 = " + frm.format(yDegerleri.get(0)));
                bDegerleri.add(yDegerleri.get(0));
            }
        }
        double sonuc = 0.0;
        double xCarpim = 1.0;
        for (int i = 0; i < bDegerleri.size(); i++) {
            for (int j = xDegerleri.size(); j > xDegerleri.size() - i; j--) {
                xCarpim *= (soruX - xDegerleri.get(xDegerleri.size() - j));
            }
            sonuc += bDegerleri.get(i) * xCarpim;
            xCarpim = 1.0;
        }
        System.out.println("f(" + soruX + ") = " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void lagrangeInterpolasyonHesapla(ArrayList<Double> xDegerleri, ArrayList<Double> yDegerleri, double soruX, double gercekDeger) {
        double sonuc = 0.0;
        double xCarpim = 1.0, payCarpim = 1.0, paydaCarpim = 1.0;
        for (int i = 0; i < yDegerleri.size(); i++) {
            for (int j = 0; j < xDegerleri.size(); j++) {
                if (j == i) {
                    continue;
                }
                paydaCarpim *= (xDegerleri.get(i) - xDegerleri.get(j));
                payCarpim *= (soruX - xDegerleri.get(j));
            }
            System.out.println("L"+i+" (x) = "+frm.format((payCarpim / paydaCarpim)));
            sonuc += (payCarpim / paydaCarpim) * yDegerleri.get(i);
            payCarpim = 1.0;
            paydaCarpim = 1.0;
        }
        System.out.println("f(" + soruX + ") = " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void fonksiyondanDegerHesapla(double x0Degeri, double xNDegeri, double h) {
        int i = 0;
        while (true) {
            double yDegeri =Math.cos((Math.PI*x0Degeri)/(6.0));
            System.out.println("x(" + i + ") = " + frm.format(x0Degeri) + "\t y(" + i + ") = " + frm.format(yDegeri));
            i++;
            x0Degeri += h;
            if (x0Degeri > xNDegeri) {
                break;
            }
        }
    }

    void ileriNumerikTurevHesapla(ArrayList<Double> deltaYLer, double h, double gercekDeger) {
        double parantezIci = 0.0;
        for (int i = 0; i < deltaYLer.size(); i++) {
            if (i % 2 == 1) {
                parantezIci -= deltaYLer.get(i) / (i + 1);
            } else {
                parantezIci += deltaYLer.get(i) / (i + 1);
            }
        }
        double sonuc = (1.0 / h) * (parantezIci);
        System.out.println("Ileri Numerik Turev x : " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void geriNumerikTurevHesapla(ArrayList<Double> deltaYLer, double h, double gercekDeger) {
        double parantezIci = 0.0;
        for (int i = 0; i < deltaYLer.size(); i++) {
            parantezIci += deltaYLer.get(i) / (i + 1);
        }
        double sonuc = (1.0 / h) * (parantezIci);
        System.out.println("Geri Numerik Turev x : " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void merkeziNumerikTurevHesapla(double ileriNTurev, double geriNTurev, double gercekDeger) {
        double sonuc = (ileriNTurev + geriNTurev) / 2.0;
        System.out.println("Merkezi Numerik Turev x : " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void yamukKuraliHesapla(ArrayList<Double> yDegerleri, double h, double gercekDeger) {
        double parantezIci = 0.0;
        for (int i = 0; i < yDegerleri.size(); i++) {
            if (i == 0 || i == yDegerleri.size() - 1) {
                parantezIci += yDegerleri.get(i);
            } else {
                parantezIci += 2 * yDegerleri.get(i);
            }
        }
        double sonuc = (h / 2.0) * (parantezIci);
        System.out.println("Yamuk Kurali x Yaklasik Olarak : " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void simpson1Bölü3Hespla(ArrayList<Double> yDegerleri, double h, double gercekDeger) {
        double parantezIci = 0.0;
        for (int i = 0; i < yDegerleri.size(); i++) {
            if (i == 0 || i == yDegerleri.size() - 1) {
                parantezIci += yDegerleri.get(i);
            } else {
                if (i % 2 == 0) {
                    parantezIci += 2 * yDegerleri.get(i);
                } else {
                    parantezIci += 4 * yDegerleri.get(i);
                }
            }
        }
        double sonuc = (h / 3.0) * (parantezIci);
        System.out.println("Simpson 1 bölü 3 x Yaklasik Olarak : " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void simpson3Bölü8Hespla(ArrayList<Double> yDegerleri, double h, double gercekDeger) {
        double parantezIci = 0.0;
        for (int i = 0; i < yDegerleri.size(); i++) {
            if (i == 0 || i == yDegerleri.size() - 1) {
                parantezIci += yDegerleri.get(i);
            } else {
                if (i % 3 == 0) {
                    parantezIci += 2 * yDegerleri.get(i);
                } else {
                    parantezIci += 3 * yDegerleri.get(i);
                }
            }
        }
        double sonuc = ((3 * h) / 8.0) * (parantezIci);
        System.out.println("Simpson 3 bölü 8 x Yaklasik Olarak : " + frm.format(sonuc));
        double hata = gercekDeger - sonuc;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - sonuc) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void ikinciRungeKuttaHesapla(double x0, double y0, double soruX, double h, double gercekDeger) {
        double x1, y1 = y0, k1, k2;
        double i = x0;
        while (true) {
            if (i == soruX) {
                break;
            } else {
                // k1 k2 verilen fonksiyon(türev dy/dx) a eşit değiştirmeyi unutma !
                k1 = (-2 * y0) / x0;
                x0 += h;
                y0 += h * k1;
                k2 = (-2 * y0) / x0;
                y0 = y1 + ((h / 2.0) * (k1 + k2));
                y1 = y0;
                i += h;
                System.out.println("y(" + (i) + ") = " + frm.format(y0));
            }
        }
        double hata = gercekDeger - y0;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - y0) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void ucuncuRungeKuttaHesapla(double x0, double y0, double soruX, double h, double gercekDeger) {
        double x1, y1 = y0, k1, k2, k3;
        double ilkDegerX = x0, ilkDegerY = y0;
        double i = x0;
        while (true) {
            if (i == soruX) {
                break;
            } else {
                k1 = ((2.0 * x0) - y0) / (1.0 + x0);
                x0 += h / 2.0;
                y0 += (h * k1) / 2.0;
                k2 = ((2.0 * x0) - y0) / (1.0 + x0);
                x0 = ilkDegerX;
                y0 = ilkDegerY;
                x0 += h;
                y0 = y0 - (h * k1) + (2.0 * h * k2);
                k3 = ((2.0 * x0) - y0) / (1.0 + x0);
                y0 = y1 + ((h / 6.0) * (k1 + 4.0 * k2 + k3));
                y1 = y0;
                ilkDegerX = x0;
                ilkDegerY = y0;
                i += h;
                System.out.println("y(" + (i) + ") = " + frm.format(y0));
            }
        }
        double hata = gercekDeger - y0;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - y0) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void dorduncuRungeKuttaHesapla(double x0, double y0, double soruX, double h, double gercekDeger) {
        double x1, y1 = y0, k1, k2, k3, k4;
        double ilkDegerX = x0, ilkDegerY = y0;
        double i = x0;
        while (true) {
            if (i == soruX) {
                break;
            } else {
                k1 = (1.0 + (x0 * y0)) / y0;
                x0 += h / 2.0;
                y0 += (h * k1) / 2.0;
                k2 = (1.0 + (x0 * y0)) / y0;
                x0 = ilkDegerX;
                y0 = ilkDegerY;
                x0 += h / 2.0;
                y0 += (h * k2) / 2.0;
                k3 = (1.0 + (x0 * y0)) / y0;
                x0 = ilkDegerX;
                y0 = ilkDegerY;
                x0 += h;
                y0 += (h * k3);
                k4 = (1.0 + (x0 * y0)) / y0;
                y0 = y1 + ((h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4));
                y1 = y0;
                ilkDegerX = x0;
                ilkDegerY = y0;
                i += h;
                System.out.println("y(" + (i) + ") = " + frm.format(y0));
            }
        }
        double hata = gercekDeger - y0;
        System.out.println("e = " + frm.format(hata));
        double bagilHata = ((gercekDeger - y0) / gercekDeger) * 100;
        System.out.println("Bagil Hata = %" + frm.format(bagilHata));
    }

    void dogrusalRegresyonHesapla(ArrayList<Double> xDegerleri, ArrayList<Double> yDegerleri) {
        double n = (double) xDegerleri.size();
        double EX = 0.0, EX2 = 0.0, EY = 0.0, EXY = 0.0;
        double a0, a1;
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX += xDegerleri.get(i);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX2 += Math.pow(xDegerleri.get(i), 2.0);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EY += yDegerleri.get(i);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EXY += xDegerleri.get(i) * yDegerleri.get(i);
        }
        double detA = (EX2 * n) - (EX * EX);
        double detA0 = (EX2 * EY) - (EXY * EX);
        double detA1 = (EXY * n) - (EX * EY);
        System.out.println("+det A :"+frm.format(detA)+"\tdet A0 :"+frm.format(detA0)+"\tdet A1 :"+frm.format(detA1));
        a0 = detA0 / detA;
        a1 = detA1 / detA;
        System.out.println("a0 : " + frm.format(a0) + " a1 : " + frm.format(a1));
    }

    void ikinciDereceRegresyonHesapla(ArrayList<Double> xDegerleri, ArrayList<Double> yDegerleri) {
        double n = (double) xDegerleri.size();
        double EX = 0.0, EX2 = 0.0, EX3 = 0.0, EX4 = 0.0, EY = 0.0, EXY = 0.0, EX2Y = 0.0;
        double a0, a1, a2;
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX += xDegerleri.get(i);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX2 += Math.pow(xDegerleri.get(i), 2.0);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX3 += Math.pow(xDegerleri.get(i), 3.0);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX4 += Math.pow(xDegerleri.get(i), 4.0);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EY += yDegerleri.get(i);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EXY += xDegerleri.get(i) * yDegerleri.get(i);
        }
        for (int i = 0; i < xDegerleri.size(); i++) {
            EX2Y += Math.pow(xDegerleri.get(i), 2.0) * yDegerleri.get(i);
        }
        ArrayList<Double> xDegerForDet = new ArrayList<>();
        xDegerForDet.add(n);
        xDegerForDet.add(EX);
        xDegerForDet.add(EX2);
        xDegerForDet.add(EX);
        xDegerForDet.add(EX2);
        xDegerForDet.add(EX3);
        xDegerForDet.add(EX2);
        xDegerForDet.add(EX3);
        xDegerForDet.add(EX4);
        ArrayList<Double> yDegerForDet = new ArrayList<>();
        yDegerForDet.add(EY);
        yDegerForDet.add(EXY);
        yDegerForDet.add(EX2Y);
        double detA = determinantHesaplama(xDegerForDet);
        double detA0 = v1MatrisHesapla(xDegerForDet, yDegerForDet);
        double detA1 = v2MatrisHesapla(xDegerForDet, yDegerForDet);
        double detA2 = v3MatrisHesapla(xDegerForDet, yDegerForDet);
        System.out.println("+det A :"+frm.format(detA)+"\tdet A0 :"+frm.format(detA0)+"\tdet A1 :"+frm.format(detA1)+"\tdet A2 :"+frm.format(detA2));
        a0 = detA0 / detA;
        a1 = detA1 / detA;
        a2 = detA2 / detA;

        System.out.println("a0 : " + frm.format(a0) + " a1 : " + frm.format(a1) + " a2 : " + frm.format(a2));
    }

    double determinantHesaplama(ArrayList<Double> matris) {
        double minorA11, minorA12, minorA13;
        minorA11 = ((double) matris.get(4) * (double) matris.get(8)) - ((double) matris.get(5) * (double) matris.get(7));
        minorA12 = ((double) matris.get(3) * (double) matris.get(8)) - ((double) matris.get(5) * (double) matris.get(6));
        minorA13 = ((double) matris.get(3) * (double) matris.get(7)) - ((double) matris.get(4) * (double) matris.get(6));
        minorA12 = -1 * minorA12;
        double kofA11, kofA12, kofA13;
        kofA11 = minorA11 * (double) matris.get(0);
        kofA12 = minorA12 * (double) matris.get(1);
        kofA13 = minorA13 * (double) matris.get(2);
        double determinant = kofA11 + kofA12 + kofA13;
        return determinant;

    }

    double v1MatrisHesapla(ArrayList<Double> v1Matris, ArrayList<Double> cozumMatrisi) {
        ArrayList<Double> v1Kopya = (ArrayList<Double>) v1Matris.clone();
        int j = 0;
        for (int i = 0; i < 9; i += 3) {
            v1Kopya.set(i, cozumMatrisi.get(j));
            j++;
        }
        return determinantHesaplama(v1Kopya);
    }

    double v2MatrisHesapla(ArrayList<Double> v2Matris, ArrayList<Double> cozumMatrisi) {
        ArrayList<Double> v2Kopya = (ArrayList<Double>) v2Matris.clone();
        int j = 0;
        for (int i = 0; i < 9; i += 3) {
            v2Kopya.set(i + 1, cozumMatrisi.get(j));
            j++;
        }
        return determinantHesaplama(v2Kopya);
    }

    double v3MatrisHesapla(ArrayList<Double> v3Matris, ArrayList<Double> cozumMatrisi) {
        ArrayList<Double> v3Kopya = (ArrayList<Double>) v3Matris.clone();
        int j = 0;
        for (int i = 0; i < 9; i += 3) {
            v3Kopya.set(i + 2, cozumMatrisi.get(j));
            j++;
        }
        return determinantHesaplama(v3Kopya);
    }
}
