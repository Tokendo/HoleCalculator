using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.Linq;
using System.Numerics;

List<double> patternHoles = new List<double> { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
    100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195,
    200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295,
    300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, };

List<double> freeHoles = new List<double>();

List<List<double>> DrillConsigne = new List<List<double>>();
double[,] drillBit = new double[1000, 3];
const double masseVolumique = 8300;
const double rayon = 5; // mm
const double degToRad = Math.PI / 180, radToDeg = 180 / Math.PI;
const long TicksPerMicrosecond = 10;

Stopwatch sw = Stopwatch.StartNew();

#region setDrillBit
double totalVolRemoved = 0, roundIntegral = 2 * Math.PI;
for (int i = 0; i < 1000; i++)
{
    drillBit[i, 0] = i; //depth in mu
    drillBit[i, 1] = (double)i / 2;//diametre in mu
    totalVolRemoved += roundIntegral * Math.Pow((drillBit[i, 1] / 2), 2);
    drillBit[i, 2] = totalVolRemoved; // vol in mu3
}
#endregion

double balourd = 50;
double angleBalourd = 180;

double balourdCible = 0;
double angleCible = 0;
double inertiaToRemove = 0;

double maxHoleDepth = 200;//in mu
double maxAngleSpan = 120, halfmaxAngleSpan = maxAngleSpan / 2;
double maxMassPerHole = getMassFromHoleDepth(maxHoleDepth);

# region methodes
double getMassFromHoleDepth(double hole)
{
    for (int i = 0; i < 1000; i++)
    {
        if (drillBit[i, 0] >= hole)
        {
            return drillBit[i, 2] * masseVolumique * Math.Pow(10, -9);//Microgramme per hole
        }
    }
    return -1;
}
double getHoleDepthFromMass(double mass)
{
    mass = mass / masseVolumique * Math.Pow(10, 9);
    for (int i = 0; i < 1000; i++)
    {
        if (drillBit[i, 2] >= mass)
        {
            return drillBit[i, 0];//micrometer
        }
    }
    return -1;
}
void getUsinageBrutAmplitude(double inputBalourd, double inputAngleDeg, double balourdCible, double angleCible, out double amplitudeBrut, out double angleBrutDeg)
{
    double balourdX = 0, balourdY = 0, cibleX = 0, cibleY = 0, usinageX = 0, usinageY = 0;

    balourdX = inputBalourd * Math.Cos(inputAngleDeg * degToRad);
    balourdY = inputBalourd * Math.Sin(inputAngleDeg * degToRad);

    cibleX = balourdCible * Math.Cos(angleCible * degToRad);
    cibleY = balourdCible * Math.Sin(angleCible * degToRad);

    usinageX = balourdX - cibleX;
    usinageY = balourdY - cibleY;

    amplitudeBrut = Math.Round(Math.Sqrt(usinageX * usinageX + usinageY * usinageY), 3);
    angleBrutDeg = Math.Round(radToDeg * Math.Atan2(usinageY, usinageX), 3);
    while (angleBrutDeg < 0)
    { angleBrutDeg += 360; }
}

void getCorrectedUsinage(double inputBalourd, double inputAngleDeg, out double amplitudeCorected, out double angleCorectedDeg)
{
    double correctorAngle = 0, correctorUsinage = 1;

    amplitudeCorected = Math.Round(inputBalourd * correctorUsinage, 3);
    angleCorectedDeg = Math.Round(inputAngleDeg + correctorAngle, 3);

    while (angleCorectedDeg < 0)
    { angleCorectedDeg += 360; }

    while (angleCorectedDeg >= 360)
    { angleCorectedDeg -= 360; }
}

double GetInertiaRemoved(List<List<double>> holesToDrill)
{
    double inertiaRemoved = 0;
    double rayonCm = rayon / 10;
    foreach (List<double> hole in holesToDrill)
    {
        inertiaRemoved += rayonCm * rayonCm * getMassFromHoleDepth(hole[1]);
    }
    return inertiaRemoved;
}

bool computeDrilling(double amplitudeUsi, double angleUsi, double inertia, out List<List<double>> holesToDrill)
{
    Console.WriteLine("Amplitude:" + amplitudeUsi.ToString() + " angle:" + angleUsi.ToString());
    List<List<double>> leftFreeHoles = new List<List<double>>(); //row { angle, diffAngle }
    List<List<double>> rightFreeHoles = new List<List<double>>(); //row { angle, diffAngle }
    holesToDrill = new List<List<double>>();

    double diffAngle = 0;
    // set possible holes on each side
    foreach (double angle in freeHoles)
    {
        diffAngle = angleUsi - angle;

        if (diffAngle < -180)
            diffAngle += 360;

        if (diffAngle > 180)
            diffAngle -= 360;

        if (Math.Abs(diffAngle) < halfmaxAngleSpan)
        {
            List<double> hole = new List<double>() { angle, diffAngle };
            if (diffAngle <= 0)
            {
                hole[1] = -hole[1];
                leftFreeHoles.Add(hole);
            }
            else
            {
                rightFreeHoles.Add(hole);
            }
        }
    }
    if (leftFreeHoles.Count == 0 || rightFreeHoles.Count == 0)
        return false;

    //Low angle difference with balourd are more efficient
    leftFreeHoles.Sort((a, b) => a[1].CompareTo(b[1]));
    rightFreeHoles.Sort((a, b) => a[1].CompareTo(b[1]));

    List<double> vRestBaloursXY = new List<double>() { amplitudeUsi * Math.Cos(angleUsi * degToRad), amplitudeUsi * Math.Sin(angleUsi * degToRad), amplitudeUsi };

    List<double> holesToDrillLeft = new List<double>(); //row angle
    holesToDrillLeft.Add(leftFreeHoles[0][0]);
    leftFreeHoles.RemoveAt(0);

    List<double> holesToDrillRight = new List<double>(); //row {angle, depth}
    holesToDrillRight.Add(rightFreeHoles[0][0]);
    rightFreeHoles.RemoveAt(0);

    bool addHoleOnLeft = true;

    while (!updateHolesDepth(amplitudeUsi, angleUsi, holesToDrillLeft, holesToDrillRight, out holesToDrill, out addHoleOnLeft))
    {

        if (addHoleOnLeft)
        {
            if (leftFreeHoles.Count != 0)
            {
                holesToDrillLeft.Add(leftFreeHoles[0][0]);
                leftFreeHoles.RemoveAt(0);
            }
            else
            { return false; }
        }
        else
        {
            if (rightFreeHoles.Count != 0)
            {
                holesToDrillRight.Add(rightFreeHoles[0][0]);
                rightFreeHoles.RemoveAt(0);
            }
            else
            { return false; }
        }
    }
    //possible to balance

    if (inertiaToRemove > 0)
    {

        double inertiaDrilled = GetInertiaRemoved(holesToDrill);
        if (inertiaDrilled >= inertiaToRemove)
        {
            return true;
        }
        else // add hole to remove inertia
        {
            List<List<double>> oppositLeftFreeHoles = new List<List<double>>(); //row { angle, diffAngle }
            List<List<double>> oppositRightFreeHoles = new List<List<double>>(); //row { angle, diffAngle }
            foreach (double angle in freeHoles)
            {
                diffAngle = angleUsi - angle;

                if (diffAngle < -180)
                    diffAngle += 360;

                if (Math.Abs(diffAngle) > halfmaxAngleSpan)
                {
                    List<double> hole = new List<double>() { angle, diffAngle };
                    if (diffAngle <= 0)
                    {
                        hole[1] = -hole[1];
                        oppositLeftFreeHoles.Add(hole);
                    }
                    else
                    {
                        oppositRightFreeHoles.Add(hole);
                    }
                }
            }
            if (oppositLeftFreeHoles.Count == 0 && oppositRightFreeHoles.Count == 0)
                return false;
            oppositLeftFreeHoles.Sort((a, b) => b[1].CompareTo(a[1]));
            oppositRightFreeHoles.Sort((a, b) => b[1].CompareTo(a[1]));


            List<double> holesToDrillOpposit = new List<double>(); //row angle
            if (oppositLeftFreeHoles[0][1] > oppositRightFreeHoles[0][1]) // get the oppositer angle
            {
                holesToDrillOpposit.Add(oppositLeftFreeHoles[0][0]);
                oppositLeftFreeHoles.RemoveAt(0);
            }
            else
            {
                holesToDrillOpposit.Add(oppositRightFreeHoles[0][0]);
                oppositRightFreeHoles.RemoveAt(0);
            }


            int state = 0;
            while (true)
            {
                state = updateHolesDepthWithInertia(amplitudeUsi, angleUsi, inertiaToRemove, holesToDrillLeft, holesToDrillRight, holesToDrillOpposit, out holesToDrill, out addHoleOnLeft);

                switch (state)
                {
                    case 0: // balanced
                        return true;
                        break;

                    case 1: // add hole imbalance side
                        if (addHoleOnLeft)
                        {
                            if (leftFreeHoles.Count != 0)
                            {
                                holesToDrillLeft.Add(leftFreeHoles[0][0]);
                                leftFreeHoles.RemoveAt(0);
                            }
                            else
                            { return false; }
                        }
                        else
                        {
                            if (rightFreeHoles.Count != 0)
                            {
                                holesToDrillRight.Add(rightFreeHoles[0][0]);
                                rightFreeHoles.RemoveAt(0);
                            }
                            else
                            { return false; }
                        }
                        break;

                    case 2: // add hole opposit side
                        if (addHoleOnLeft)
                        {
                            if (oppositLeftFreeHoles.Count != 0)
                            {
                                holesToDrillOpposit.Add(oppositLeftFreeHoles[0][0]);
                                oppositLeftFreeHoles.RemoveAt(0);
                            }
                            else
                            { return false; }
                        }
                        else
                        {
                            if (oppositRightFreeHoles.Count != 0)
                            {
                                holesToDrillOpposit.Add(oppositRightFreeHoles[0][0]);
                                oppositRightFreeHoles.RemoveAt(0);
                            }
                            else
                            { return false; }
                        }
                        break;

                    case 3: // balancing not possible
                        return false;
                }
            }
        }

    }
    else
        return true;
}
bool updateHolesDepth(double amplitudeUsi, double angleUsi, List<double> holesToDrillLeft, List<double> holesToDrillRight, out List<List<double>> holesToDrillUpdated, out bool addHoleLeft)
{
    holesToDrillUpdated = new List<List<double>>();
    addHoleLeft = true;
    double usinageCibleX = amplitudeUsi * Math.Cos(angleUsi * degToRad);
    double usinageCibleY = amplitudeUsi * Math.Sin(angleUsi * degToRad);
    double rayonCm = rayon / 10;
    double usinageFullHolesX = 0, usinageFullHolesY = 0;


    //Firsts hole are drilled to max, we need to remove mass
    if (holesToDrillLeft.Count > 1)
    {
        for (int i = 0; i < holesToDrillLeft.Count - 1; i++)
        {
            usinageFullHolesX += maxMassPerHole * rayonCm * Math.Cos(holesToDrillLeft[i] * degToRad);
            usinageFullHolesY += maxMassPerHole * rayonCm * Math.Sin(holesToDrillLeft[i] * degToRad);
            holesToDrillUpdated.Add(new List<double> { holesToDrillLeft[i], maxHoleDepth });
        }
    }

    if (holesToDrillRight.Count > 1)
    {
        for (int i = 0; i < holesToDrillRight.Count - 1; i++)
        {
            usinageFullHolesX += maxMassPerHole * rayonCm * Math.Cos(holesToDrillRight[i] * degToRad);
            usinageFullHolesY += maxMassPerHole * rayonCm * Math.Sin(holesToDrillRight[i] * degToRad);
            holesToDrillUpdated.Add(new List<double> { holesToDrillRight[i], maxHoleDepth });
        }
    }

    usinageCibleX = usinageCibleX - usinageFullHolesX;
    usinageCibleY = usinageCibleY - usinageFullHolesY;

    //The last 2 holes balance the drilling to get in target
    // cibleX= gaucheX*a + droitX*b
    // cibleY= gaucheY*a + droitY*b
    // a = masse perçage gauche
    // b = masse perçage droit
    // a et b < maxHoleDepth pour ok
    double a = 0, b = 0;
    double gaucheX = rayonCm * Math.Cos(holesToDrillLeft.Last() * degToRad);
    double gaucheY = rayonCm * Math.Sin(holesToDrillLeft.Last() * degToRad);

    double droitX = rayonCm * Math.Cos(holesToDrillRight.Last() * degToRad);
    double droitY = rayonCm * Math.Sin(holesToDrillRight.Last() * degToRad);

    b = (usinageCibleY - (gaucheY / gaucheX) * usinageCibleX) / ((-gaucheY / gaucheX) * droitX + droitY);
    a = (usinageCibleX - (droitX * b)) / gaucheX;

    if (a >= 0 && a <= maxMassPerHole && b >= 0 && b <= maxMassPerHole)
    {
        //convert mass to hole depth
        a = getHoleDepthFromMass(a);
        b = getHoleDepthFromMass(b);
        holesToDrillUpdated.Add(new List<double> { holesToDrillLeft.Last(), a });
        holesToDrillUpdated.Add(new List<double> { holesToDrillRight.Last(), b });
        Console.WriteLine("Solution trouvée: prof a=" + a.ToString("F3") + " prof b=" + b.ToString("F3"));
        return true;// equillibrage ok
    }
    else
    {
        addHoleLeft = a >= b;
        return false;//add equillibrage hole
    }
}

int updateHolesDepthWithInertia(double amplitudeUsi, double angleUsi, double inertia,
    List<double> holesToDrillLeft, List<double> holesToDrillRight,
    List<double> holesToDrillOpposit,
    out List<List<double>> holesToDrillUpdated, out bool addHoleLeft)
{
    holesToDrillUpdated = new List<List<double>>();
    addHoleLeft = true;
    double usinageCibleX = amplitudeUsi * Math.Cos(angleUsi * degToRad);
    double usinageCibleY = amplitudeUsi * Math.Sin(angleUsi * degToRad);
    double rayonCm = rayon / 10;
    double usinageFullHolesX = 0, usinageFullHolesY = 0;
    double inertiaRemoved = 0;


    //Firsts hole are drilled to max, we need to remove mass
    if (holesToDrillLeft.Count > 1)
    {
        for (int i = 0; i < holesToDrillLeft.Count - 1; i++)
        {
            usinageFullHolesX += maxMassPerHole * rayonCm * Math.Cos(holesToDrillLeft[i] * degToRad);
            usinageFullHolesY += maxMassPerHole * rayonCm * Math.Sin(holesToDrillLeft[i] * degToRad);
            holesToDrillUpdated.Add(new List<double> { holesToDrillLeft[i], maxHoleDepth });
            inertiaRemoved += rayonCm * rayonCm * maxMassPerHole;
        }
    }

    if (holesToDrillRight.Count > 1)
    {
        for (int i = 0; i < holesToDrillRight.Count - 1; i++)
        {
            usinageFullHolesX += maxMassPerHole * rayonCm * Math.Cos(holesToDrillRight[i] * degToRad);
            usinageFullHolesY += maxMassPerHole * rayonCm * Math.Sin(holesToDrillRight[i] * degToRad);
            holesToDrillUpdated.Add(new List<double> { holesToDrillRight[i], maxHoleDepth });
            inertiaRemoved += rayonCm * rayonCm * maxMassPerHole;
        }
    }

    if (holesToDrillOpposit.Count > 1)
    {
        for (int i = 0; i < holesToDrillOpposit.Count - 1; i++)
        {
            usinageFullHolesX += maxMassPerHole * rayonCm * Math.Cos(holesToDrillOpposit[i] * degToRad);
            usinageFullHolesY += maxMassPerHole * rayonCm * Math.Sin(holesToDrillOpposit[i] * degToRad);
            holesToDrillUpdated.Add(new List<double> { holesToDrillOpposit[i], maxHoleDepth });
            inertiaRemoved += rayonCm * rayonCm * maxMassPerHole;
        }
    }

    usinageCibleX = usinageCibleX - usinageFullHolesX;
    usinageCibleY = usinageCibleY - usinageFullHolesY;

    // les 3 derniers trous doivent équilibrer et arriver à la cible d'inertie
    //cibleX=gaucheX*a+droitX*b+opposeX*c
    //cibleY=gaucheY*a+droitY*b+opposeY*c
    //inertia=rayon^2*(a+b+c) ==> inertia/rayon^2=a+b+c

    double a = 0, b = 0, c = 0;
    double gaucheX = rayonCm * Math.Cos(holesToDrillLeft.Last() * degToRad);
    double gaucheY = rayonCm * Math.Sin(holesToDrillLeft.Last() * degToRad);

    double droitX = rayonCm * Math.Cos(holesToDrillRight.Last() * degToRad);
    double droitY = rayonCm * Math.Sin(holesToDrillRight.Last() * degToRad);

    double oppX = rayonCm * Math.Cos(holesToDrillOpposit.Last() * degToRad);
    double oppY = rayonCm * Math.Sin(holesToDrillOpposit.Last() * degToRad);

    double lastInertia = Convert.ToDouble(inertia) - inertiaRemoved;

    MathNet.Numerics.LinearAlgebra.Matrix<double> Amat = Matrix<double>.Build.DenseOfArray(new double[,]
    {
    { gaucheX, droitX, oppX },
    { gaucheY, droitY, oppY },
    { 1, 1, 1 }
    });
    MathNet.Numerics.LinearAlgebra.Vector<double> bv = MathNet.Numerics.LinearAlgebra.Vector<double>.Build.Dense(new double[] { usinageCibleX, usinageCibleY, lastInertia / (rayonCm * rayonCm) });
    var x = Amat.Solve(bv);
    a = x[0]; b = x[1]; c = x[2];

    Console.WriteLine("masse a=" + a.ToString("F3") + " b=" + b.ToString("F3") + " c=" + c.ToString("F3"));

    if (a >= 0 && a <= maxMassPerHole && b >= 0 && b <= maxMassPerHole && c >= 0 && c <= maxMassPerHole)
    {
        //convert mass to hole depth
        a = getHoleDepthFromMass(a);
        b = getHoleDepthFromMass(b);
        c = getHoleDepthFromMass(c);
        holesToDrillUpdated.Add(new List<double> { holesToDrillLeft.Last(), a });
        holesToDrillUpdated.Add(new List<double> { holesToDrillRight.Last(), b });
        holesToDrillUpdated.Add(new List<double> { holesToDrillOpposit.Last(), c });
        Console.WriteLine("Solution trouvée, prof a=" + a.ToString("F3") + " prof b=" + b.ToString("F3") + " prof c=" + c.ToString("F3"));
        return 0;// equillibrage ok
    }
    else if (c > a && c > b) //need more opposit
    {
        addHoleLeft = a >= b;
        return 2;//add opposit hole
    }
    else
    {
        addHoleLeft = a >= b;
        return 1;//add opposit hole
    }
}
void drawUsinage(double angleBalourd, List<List<double>> consigneUsinage, List<double> holesDisp)
{
    int bmpSize = 1000;
    int center = bmpSize / 2;
    int radius = bmpSize / 3;

    int x = 0, y = 0;
    var bmp = new Bitmap(bmpSize, bmpSize);
    Graphics g = Graphics.FromImage(bmp);
    Pen pen = new Pen(Color.Gold);
    Point startPoint = new Point(0, 0);
    Point endPoint = new Point(0, 0);
    Point balCenter = new Point(center, center);
    Rectangle rect = new Rectangle(center - radius, center - radius, radius * 2, radius * 2);
    g.DrawEllipse(pen, rect);

    pen = new Pen(Color.White);

    x = Convert.ToInt32(1.0 * center + 1.2 * radius * Math.Cos(0));
    y = Convert.ToInt32(1.0 * center - 1.2 * radius * Math.Sin(0));
    endPoint = new Point(x, y);
    g.DrawLine(pen, balCenter, endPoint);

    // draw freeholes
    pen = new Pen(Color.Green);
    for (int i = 0; i < holesDisp.Count; i++)
    {
        x = Convert.ToInt32(1.0 * center + (double)radius * Math.Cos(holesDisp[i] * degToRad));
        y = Convert.ToInt32(1.0 * center - (double)radius * Math.Sin(holesDisp[i] * degToRad));
        endPoint = new Point(x, y);
        g.DrawLine(pen, balCenter, endPoint);
    }

    // draw drilled holes
    pen = new Pen(Color.Blue);
    for (int i = 0; i < consigneUsinage.Count; i++)
    {
        double amplitud = 0;
        amplitud = consigneUsinage[i][1] / maxHoleDepth;

        // Console.WriteLine("Angle: " + consigneUsinage[i][0].ToString() + " Amplitude:" + amplitud.ToString() + " Profondeur:" + consigneUsinage[i][1].ToString());

        x = Convert.ToInt32(1.0 * center + (amplitud * radius) * Math.Cos(consigneUsinage[i][0] * degToRad));
        y = Convert.ToInt32(1.0 * center - (amplitud * radius) * Math.Sin(consigneUsinage[i][0] * degToRad));
        endPoint = new Point(x, y);
        g.DrawLine(pen, balCenter, endPoint);
    }

    //draw imbalance
    pen = new Pen(Color.Violet);
    x = Convert.ToInt32(1.0 * center + (10.2 + radius) * Math.Cos(angleBalourd * degToRad));
    y = Convert.ToInt32(1.0 * center - (10.2 + radius) * Math.Sin(angleBalourd * degToRad));
    endPoint = new Point(x, y);
    g.DrawLine(pen, balCenter, endPoint);

    // draw zero


    g.Dispose();
    bmp.Save("Balancier balourd" + balourd.ToString("F3") + " , angle" + angleBalourd.ToString("F3") + ".png", ImageFormat.Png);
}

#endregion
Random rand = new Random();
for (int i = 0; i < 200; i++)
{
    freeHoles.Clear();
    foreach (double hole in patternHoles)
    {
        if (rand.NextDouble() > 0.2)
            freeHoles.Add(hole);
    }
    sw.Restart();
    balourd = 100.0 * rand.NextDouble();
    angleBalourd = 360 * rand.NextDouble();

    double amplitudeUsinageBrut = 0, angleUsinageBrut = 0;

    getUsinageBrutAmplitude(balourd, angleBalourd, balourdCible, angleCible, out amplitudeUsinageBrut, out angleUsinageBrut);

    double amplitudeUsinageReel = 0, angleUsinageReel = 0;

    getCorrectedUsinage(amplitudeUsinageBrut, angleUsinageBrut, out amplitudeUsinageReel, out angleUsinageReel);

    DrillConsigne.Clear();

    if (computeDrilling(amplitudeUsinageReel, angleUsinageReel, 0, out DrillConsigne))
        //drawUsinage(angleBalourd, DrillConsigne, freeHoles);
        ;
    else
        Console.WriteLine("pas possible");
    sw.Stop();
    long microseconds = sw.ElapsedTicks / TicksPerMicrosecond;
    Console.WriteLine("duree de traitement:" + microseconds.ToString() + "microsec");
}
