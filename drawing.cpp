#ifdef WIN32
#include <windows.h>
#endif

#if defined (__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
//#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#else //linux
#include <GL/gl.h>
#include <GL/glut.h>
#endif

//other includes
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream> //need for files
#include <string>
#include <sstream>
#include <cctype> //need for isnum()
#include <cmath> //cos and sin and roundf
using namespace std;

//Global Variables!
//the number of pixels in the grid

int grid_width;
int grid_height;

int view_width = 100;
int view_height = 100;
int view_wmin = 0;
int view_hmin = 0;

//the size of pixels sets the inital window height and width
//don't make the pixels too large or the screen size will be larger than
//your display size
double pixel_size;

/*Window information*/
int win_height;
int win_width;

//global variable flags
int selectedP = -1;
int rotation = 0;
int transl_ar[2] = {0, 0};
int translx = 1;
int transly = 1;
double scale = 0.0;
double userscale = 0.2;
int exitflag = 0;
int lineAlgFlag = 1; //1 for DDA, -1 for bresenham
int angle = 45;
int showcolor = 1;

struct Polygon
{
public:
	int vertices;
	double* X;
	double* Y;
	Polygon(){} //default constructor
	Polygon(int v);
	Polygon(Polygon & rhs);
	void displayPoly(int** buffer);
	void rotateRight();
	void rotateLeft();
	void translate();
	void scalefunc();
	void scanline(int** buffer);
};



void check();
void idle();
void draw_pix(int x, int y);
void display();
void init();
void reshape(int width, int height);
void slopeZero(int x0, int y0, int x1, int y1, int** buffer);
void slopeInf(int x0, int y0, int x1, int y1, int** buffer);
void lineDDA(int x0, int y0, int x1, int y1, int** buffer);
void lineBres(int x0, int y0, int x1, int y1, int** buffer);
void bBigSlope(int x0, int y0, int x1, int y1, int** buffer, int inc);
void bSmallSlope(int x0, int y0, int x1, int y1, int** buffer, int inc);
void read_file(string* data);
void read_data(string* data, int &counter, Polygon* arrPoly, int i);
void initbuffer(int** buffer);
void freeobjects(Polygon* poly, int** buffer, string* data, int n);
void colorPoly(int** buffer);
void edge(int** buffer, int above, int below, int& i, int& j, int* status);
void motion(int x, int y);
void key(unsigned char ch, int x, int y);
void mouse(int button, int state, int x, int y);
void write_back(Polygon* poly, int numPoly);
void swapvalues(int& i, int& j);





Polygon::Polygon(int v)
{
	vertices = v;
	X = new double[v];
	Y = new double[v];
}

Polygon::Polygon(Polygon & rhs)
{
	vertices = rhs.vertices;
	for(int i = 0; i < vertices; i++)
	{
		X[i] = rhs.X[i];
		Y[i] = rhs.Y[i];		
	}
}

void swapvalues(int& i, int& j)
{
	int temp;
	temp = i;
	i = j;
	j = temp;
}


void bSmallSlope(int x0, int y0, int x1, int y1, int** buffer, int inc)
{

	int dx = fabs(x1 - x0), dy = fabs(y1 - y0);
	int p = 2 * dy - dx;
	// int p = dy - 2* dx;
	int twoDy = 2 * dy, twoDyMinusDx = 2 * (dy - dx);
	int x = x0, y = y0;

	if(x < view_width && x >= 0 && y < view_height && y >= 0)
		buffer[x][y] = 1; //plot first point

	while( x < x1)
	{
		x++;
		if(p < 0)
			p += twoDy;
		else //p >= 0
		{
			y += inc; //if inc is 1, pos m. inc is -1, neg m
			p += twoDyMinusDx;
		}

		if(x < grid_width && x >= 0 && y < grid_height && y >= 0)
			buffer[x][y] = 1;
	}
}

void bBigSlope(int x0, int y0, int x1, int y1, int** buffer, int inc)
{
	int dx = fabs(x1 - x0), dy = fabs(y1 - y0);
	int p = 2 * dx - dy;
	// int twoDy = 2 * dy, twoDyMinusDx = 2 * (dy - dx);
	int twoDx = 2 * dx, twoDxMinusDy = 2 * (dx - dy), twoDy = 2*dy;
	int x = x0, y = y0;

	if(x < grid_width && x >= 0 && y < grid_height && y >= 0)
		buffer[x][y] = 1; //plot first point

	while( y < y1)
	{
		y++;
		if(p < 0)
			p += twoDx;
		else  //p >= 0
		{
			x += inc; //if inc is 1, pos m. inc is -1, neg m
			p += twoDxMinusDy;
		}
		
		if(x < grid_width && x >= 0 && y < grid_height && y >= 0)
			buffer[x][y] = 1;
	}
}

void slopeInf(int x0, int y0, int x1, int y1, int** buffer)
{
	if(y0 > y1) 
		swapvalues(y0, y1);
		
	for(int y = y0; y < y1; y++)
		if(x0 < grid_width && x0 >= 0 && y < grid_height && y >= 0)
			buffer[x0][y] = 1;
}

void slopeZero(int x0, int y0, int x1, int y1, int** buffer)
{
	if(x0 > x1) 
		swapvalues(x0, x1);
		
	for(int x = x0; x < x1; x++)
		if(x < grid_width && x >= 0 && y0 < grid_height && y0 >= 0)
			buffer[x][y0] = 1;
}

void lineBres(int x0, int y0, int x1, int y1, int** buffer)
{
	float slope;    //  1, 20 -> 1, 12   vs   1, 12->1,20
	if(x1 - x0 == 0) //zero denom = infinite slope
	{
		slopeInf(x0, y0, x1, y1, buffer);
		return;
	}

	slope = (float(y1 - y0) / (x1 - x0));

	if(slope == 0)  //same y value
		slopeZero(x0, y0, x1, y1, buffer);

	else if(slope < 1 && slope > -1) //small slope (pos or neg)
	{
		if(x0 > x1)   //make it so that x1 is always the rightmost endpt
		{
			swapvalues(x0, x1);
			swapvalues(y0, y1);
		}

		if(slope > 0)
			bSmallSlope(x0, y0, x1, y1, buffer, 1);
		else //slope is negative
			bSmallSlope(x0, y0, x1, y1, buffer, -1);
	}

	else // slope >=1 or slope <= -1
	{
		if(y0 > y1)
		{
			swapvalues(x0, x1);
			swapvalues(y0, y1);
		}

		if(slope >= 1)
			bBigSlope(x0, y0, x1, y1, buffer, 1);
		else //negative slope
			bBigSlope(x0, y0, x1, y1, buffer, -1);
	}
		
}	//lineBres()

inline int round1(const float a) {return int (a + 0.5);}

void ddaSmallSlope(int x0, int y0, int x1, int y1, int** buffer)
{
	int dx = x1 - x0, dy = y1 - y0, steps, k;
	float xIncrement, yIncrement, x = x0, y = y0;

	if(fabs(dx) > fabs(dy))
		steps = fabs(dx);
	else
		steps = fabs(dy);

	xIncrement = float (dx) / float (steps);
	yIncrement = float (dy) / float (steps);

	if(round1(x) < grid_width && round1(x) >= 0 && round1(y) < grid_height && round1(y) >= 0)
			buffer[round1(x)][round1(y)] = 1;	
	for(k = 0; k < steps; k++)
	{
		x += xIncrement;
		y += yIncrement;
		if(round1(x) < grid_width && round1(x) >= 0 && round1(y) < grid_height && round1(y) >= 0)
			buffer[round1(x)][round1(y)] = 1;
	}
}

void lineDDA(int x0, int y0, int x1, int y1, int** buffer)
{
	float slope;    //  1, 20 -> 1, 12   vs   1, 12->1,20
	if(x1 - x0 == 0) //zero denom = infinite slope
	{
		slopeInf(x0, y0, x1, y1, buffer);
		return;
	}

	slope = (float(y1 - y0) / (x1 - x0));

	if(slope == 0)  //same y value
		slopeZero(x0, y0, x1, y1, buffer);

	else if(slope < 1 && slope > -1) //small slope (pos or neg)
	{
		if(x0 > x1)   //make it so that x1 is always the rightmost endpt
		{
			swapvalues(x0, x1);
			swapvalues(y0, y1);
		}

		if(slope > 0)
			ddaSmallSlope(x0, y0, x1, y1, buffer);
		else //slope is negative
			ddaSmallSlope(x0, y0, x1, y1, buffer);
	}

	else // slope >=1 or slope <= -1
	{
		if(y0 > y1)
		{
			swapvalues(x0, x1);
			swapvalues(y0, y1);
		}

		if(slope >= 1)
			ddaSmallSlope(x0, y0, x1, y1, buffer);
		else //negative slope
			ddaSmallSlope(x0, y0, x1, y1, buffer);
			//bBigSlope(x0, y0, x1, y1, buffer, -1);
	}
}


void Polygon::displayPoly(int** buffer)
{

	for(int v = 0; v < vertices; v++)
	{
		int x = round(X[v]);
		int y = round(Y[v]);
		if(x < grid_width && x >= 0 && 
		   y < grid_height && y >= 0 )
			buffer[x][y] = 1; //turn on the pixel if in bound
		
		if(v > 0) //theres previous vertex, then connect w line
     		if(lineAlgFlag == 1)
     			lineDDA(x, y, round(X[v-1]), round(Y[v-1]), buffer);
     		else
     			lineBres(x, y, round(X[v-1]), round(Y[v-1]), buffer); 
     	
	}

	//finally connect an edge from last vertex to first
	if(lineAlgFlag == 1)
		lineDDA(round(X[0]), round(Y[0]), round(X[vertices-1]), round(Y[vertices-1]), buffer); 
	else
		lineBres(round(X[0]), round(Y[0]), round(X[vertices-1]), round(Y[vertices-1]), buffer); 

	for(int i = view_wmin; i < view_width; i++)    // move this after color
		for(int j = view_hmin; j < view_height; j++)
			if(buffer[i][j] == 1)
				draw_pix(i, j);

} // displayPoly()


void colorPoly(int** buffer)
{
	int status = -1;
	int above = 0, below = 0, right = 0;

	for(int j = 0; j < grid_height; j++)
		for(int i = 0; i < grid_width; i++)
		{
			if(status == 1 && buffer[i][j] == 0) //reset flags if currently "in"?
			{
				below = 0;
				above = 0;
				right = 0;
				buffer[i][j] = 2; //turn on pixel (but with a 2...)
				continue;
			}

			if(status == -1 && buffer[i][j] == 0)
			{
				below = 0;
				above = 0;
				right = 0;
				continue;
			}

			if(buffer[i][j] == 1) //change status (in-1, out-(-1))
			{
				//find the bounds
				int start, end;
				if(i > 0 )
					start = i-1;
				else //i == 0
					start = i;

				if(i < grid_width-1)
					end = i+1;
				else //i == grid_width
					end = i;
				
				//for(int c = i-1; c <= i+1; c++)
				for(int c = start; c <= end; c++) //below test
					if(j-1 >= 0)
						if (buffer[c][j-1] == 1)
						{	
							below = 1;
							break;
						}
		
				for(int c = start; c <= end; c++) //above test
					if(j+1 <= grid_height-1)
						if (buffer[c][j+1] == 1)				
						{	
							above = 1;
							break;
						}		

				if(i+1 < grid_width && buffer[i+1][j] == 1)
					right = 1;

				if(right == 1) //deal with the edge
				{
					edge(buffer, above, below, i, j, &status);
					continue;
				}

				else //no adjacent pixel, check if extremity
				{
					if(above == 1 && below == 1) //normal
					{
						status *= -1;   //FIXME of course, if is more complicated
						continue;
					}

					else //above is 1 only or below is 1 only
						continue;
				}
				
			}
		}

	//flush out buffer
	for(int i = view_wmin; i < view_width; i++)    // move this after color
		for(int j = view_hmin; j < view_height; j++)
			if(buffer[i][j] == 2)
				draw_pix(i, j);

	//delete[] minmax;
}	//colorPoly()

void edge(int** buffer, int above, int below, int& i, 
				   int& j, int* status)
{
	//right flag was marked, so keep going in buffer until end,
	//then determine the end pixel's flags (above2, below2)
	int above2 = 0, below2 = 0;
	int p;
	for(p = i; p < (grid_width - 1) && buffer[p][j] == 1; p++);

	p--; //remember that p will be 1 higher than last pixel of edge
		//and after this edge func finishes, "continue" makes i++
		//which means we will skip a pixel :O
	
	int start, end;
	if(p > 0 )
		start = p-1;
	else //i == 0
		start = p;

	if(p < grid_width-1)
		end = p+1;
	else //i == grid_width
		end = p;

	for(int c = start; c <= end; c++)
		if (j-1 >= 0 && buffer[c][j-1] == 1)
		{	
			below2 = 1;
			break;
		}
				
	for(int c = start; c <= end; c++)
		if (j+1 <= grid_height && buffer[c][j+1] == 1)				
		{	
			above2 = 1;
			break;
		}	

	//update i...
	i = p;  

	//U shape
	if (above2 == 1 && above == 1 && below == 0 && below2 == 0
		|| below2 == 1 && below == 1 && above == 0 && above2 == 0)
	{
		return; //status stays the same
	}
	else if(above == 0 && below == 0) //left part edge not connected
		*status = -1;
	 

	else if (above == 1 && below == 1) //original end has both flags
		if(below2 == 1 && above2 == 1) //& end2 not connected (concave)
			*status = -1;
		else if (below2 == 1)
			*status *= -1;
		else
			*status = -1;
	else if (above2 == 1 && below2 == 1)
		*status = 1;
	else //(above2 == 1 && below == 1 || below2 == 1 && above == 1)
	{
		*status *= -1;
		
	}

	return;
} //edge()

void Polygon::translate()
{
	for(int i = 0; i < vertices; i++)
	{
		X[i] += transl_ar[0];
		Y[i] += transl_ar[1];
	}

	transl_ar[0] = 0;
	transl_ar[1] = 0;
}


// the last shot
void Polygon::rotateLeft()
{
	//actually this part is for scaling not rotation... so FIXME

	int xcen, ycen, xsum = 0, ysum = 0;

	for(int i = 0; i < vertices; i++) //find where to transl
	{
		xsum += X[i];
		ysum += Y[i];
	}

	xcen = xsum / vertices;
	ycen = ysum / vertices;

	for(int i = 0; i < vertices; i++)
	{
		//OMG huge mistake, must save previous x[i] value,
		//otherwise new x[i] will be used to calc Y[i] -_____-
		double temp;
		temp = xcen + (X[i] - xcen) * cos(angle *  -3.14159267 / 180.f)
		 		- (Y[i] - ycen) * sin(angle * -3.14159267 / 180.f);
		
		Y[i] = ycen + (X[i] - xcen) * sin(angle * -3.14159267 / 180.f)
						+ (Y[i] - ycen) * cos(angle * -3.14159267 / 180.f);
		X[i] = temp;
	}

	rotation = 0; //reset rotation

}

void Polygon::rotateRight()
{
	//actually this part is for scaling not rotation... so FIXME

	int xcen, ycen, xsum = 0, ysum = 0;

	for(int i = 0; i < vertices; i++) //find where to transl
	{
		xsum += X[i];
		ysum += Y[i];
	}

	xcen = xsum / vertices;
	ycen = ysum / vertices;

	for(int i = 0; i < vertices; i++)
	{
		double temp;
		temp = xcen + (X[i] - xcen) * cos(angle * 3.14159267 / 180.f)
		 		- (Y[i] - ycen) * sin(angle * 3.14159267 / 180.f);
		
		Y[i] = ycen + (X[i] - xcen) * sin(angle * 3.14159267 / 180.f)
						+ (Y[i] - ycen) * cos(angle * 3.14159267 / 180.f);
		X[i] = temp;
	}

	rotation = 0; //reset rotation

}


void Polygon::scalefunc()
{
	int xmin, v;
	xmin = 0;

	for(v = 0; v < vertices; v++)
		if(X[v] < X[xmin])
			xmin = v;

	double s = 1 + scale;

	for(int i = 0; i < vertices; i++)
	{
		X[i] = ((X[i] - X[xmin]) * s) + X[xmin];
		Y[i] = ((Y[i] - Y[xmin]) * s) + Y[xmin];
	}

	scale = 0.0;

		
}

void max(double* arr, double &max, int size)
{
	max = arr[0];
	for(int i = 1; i < size; i++)
		if(arr[i] > max)
			max = arr[i];
}

void min(double* arr, double &min, int size)
{
	min = arr[0];
	for(int i = 1; i < size; i++)
		if(arr[i] < min)
			min = arr[i];
}

void Polygon::scanline(int** buffer)
{
	double** edges = new double*[vertices];
	double* xIntercepts = new double[100];
	int pos = 0;
	int numedge = 0;
	double ymin, ymax, xmin, xmax;

	min(X, xmin, vertices);
	min(Y, ymin, vertices);
	max(X, xmax, vertices);
	min(Y, ymin, vertices);

	//go through each edge, during edges going up, make changes + save
	for(int v = 0; v < vertices - 1; v++)
	{
		edges[numedge] = new double[4]; //the x and y of both ends
		if(Y[v+1] != Y[v]) //dont remmeber the horiz edges
		{
			edges[numedge][0] = X[v];
			edges[numedge][1] = Y[v];
			edges[numedge][2] = X[v+1];
			edges[numedge][3] = Y[v+1];
			numedge += 1;
		}
		
	}
	if(Y[vertices - 1] != Y[0]) //dont remmeber the horiz edges
	{
		edges[numedge] = new double[4];
		edges[numedge][0] = X[vertices-1];
		edges[numedge][1] = Y[vertices-1];
		edges[numedge][2] = X[0];
		edges[numedge][3] = Y[0];
		numedge++; //now numedge == number of total edges (0->num-1)
	}

	//fix edges if needed (shorten)
	for(int e = 0; e < numedge - 1; e++)
	{
		if(edges[e][3] > edges[e+1][3]) //y values keep going up
			edges[e][3] -= 1.0;
		if(edges[e][3] < edges[e+1][3]) //y values go down 
			edges[e+1][1] -= 1.0;
	}

	if(edges[numedge - 1][3] > edges[0][3]) //y values keep going up
		edges[numedge - 1][3] -= 1.0;
	if(edges[numedge - 1][3] < edges[0][3]) //y values go down 
		edges[0][1] -= 1.0;	

	for(int j = ymin; j < ymax; j++) //each j is scanline
	{

		//find intersections, looping over edges
		//find where y = y intersects eligile edges
		// for(int e = 0; e < numedge; e++)
		// 	if(edges[e][1] >= ymin && edges[e][3] <= ymax)
				//then calculate intersection
				//check if intersection lies on next edge too
				//add intersection once or twice...
	}

}

void display()
{
	//clears the screen
	glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
    //clears the opengl Modelview transformation matrix
	glLoadIdentity();

    int numPoly, counter = 1;
    string *data = new string[100];

    read_file(data); 

    istringstream(data[0]) >> numPoly; //change string to int
    Polygon* poly = new Polygon[numPoly];

    int** buffer = new int*[grid_width];
    initbuffer(buffer);

    for(int i = 0; i < numPoly; i++) //make numPoly # of polygons
    {
    	read_data(data, counter, poly, i);
    } //read data into the polygons array, draw edges, rasterize


    if(rotation == 1)	
    	poly[selectedP].rotateRight(); //FIXME translate by true center hahah
    if(rotation == -1)
    	poly[selectedP].rotateLeft();
    if(transl_ar[0] != 0 || transl_ar[1] != 0)
    	poly[selectedP].translate();
    if(scale != 0)
    	poly[selectedP].scalefunc();

    for(int i = 0; i < numPoly; i++)
    	poly[i].displayPoly(buffer);

    // for(int i = 0; i < numPoly; i++)
    // 	poly[i].scanline(buffer);
    if(showcolor == 1)
    	colorPoly(buffer); //color all the polynomials

	write_back(poly, numPoly);
 
    //blits the current opengl framebuffer on the screen
    glutSwapBuffers();	//makes screen black except for 1 pixel
    //checks for opengl errors
	check();

	freeobjects(poly, buffer, data, numPoly);

	if(exitflag == 1)
    	exit(0);
}// display()

void write_back(Polygon* poly, int numPoly)
{
	ofstream myfile;		
	myfile.open("test_scene");	//doing this wipes the file
	myfile.precision(10);
	myfile << numPoly << endl << endl;
	for(int i = 0; i < numPoly; i++)
	{
		myfile << poly[i].vertices <<endl;
		for(int j = 0; j < poly[i].vertices; j++)
			myfile << fixed << poly[i].X[j] << " " << fixed << poly[i].Y[j] <<endl;
		myfile << endl;
	}

	myfile.close();

}

void initbuffer(int** buffer)
{
	//keep track of turned on pix
	for(int i = 0; i < grid_width; i++)
	{
		buffer[i] = new int[grid_height];
		for(int j = 0; j < grid_height; j++)
			buffer[i][j] = 0; 			//turn off all pixels
	}
}



//gets called when the curser moves accross the scene
void motion(int x, int y)
{
    //redraw the scene after mouse movement
	glutPostRedisplay();
}

//gets called when a key is pressed on the keyboard
void key(unsigned char ch, int x, int y)
{

	if(isdigit(ch))
	{
		selectedP = ch - 48;    //choose a polygon
		cout << "Selected polygon "<<selectedP<<endl;
	}

	if( ch == 'w')
		scale = userscale;
	else if( ch == 'z')
		scale = -1 * userscale;
	else if( ch == 'a')
		rotation = -1;  //left 45 degree
	else if( ch == 's')
		rotation = 1;  //left 45 degree
	else if( ch == 'i')
		transl_ar[1] = transly;
	else if( ch == 'j')
		transl_ar[0] = -1 * translx;
	else if( ch == 'k' )	
		transl_ar[1] = -1 * transly;
	else if( ch == 'l')
		transl_ar[0] = transly;
	else if (ch == 'p')
	{
		lineAlgFlag *= -1;
		if(lineAlgFlag == -1)
			cout << "p: used Bresenham"<<endl;
		else
			cout <<"p: used DDA"<<endl;
	}
	else if(ch == 'c')
		showcolor *= -1;

	else 
		scale = 0;

	switch(ch)
	{
		case 'q':
			exitflag = 1;
			cout << "q: quit"<<endl;
		// case 'a':
		// 	rotation = -1;  //left 45 degree
		// case 's':
		// 	rotation = 1;
		// case 'i':
			// transl_ar[1] = 1;
		// case'j':
		// 	transl_ar[0] = -1;
		// case 'k':	
		// 	transl_ar[1] = -1;
		// case 'l':
		// 	transl_ar[0] = 1;

		default:
            //prints out which key the user hit
            // printf("User hit the \"%c\" key\n",ch);
			break;

	}
    //redraw the scene after keyboard input
	glutPostRedisplay();
}


//gets called when a mouse button is pressed
void mouse(int button, int state, int x, int y)
{
    //print the pixel location, and the grid location
    printf ("MOUSE AT PIXEL: %d %d, GRID: %d %d\n",x,y,(int)(x/pixel_size),(int)(y/pixel_size));
	switch(button)
	{
		case GLUT_LEFT_BUTTON: //left button
            printf("LEFT ");
            break;
		case GLUT_RIGHT_BUTTON: //right button
            printf("RIGHT ");
		default:
            printf("UNKNOWN "); //any other mouse button
			break;
	}
    if(state !=GLUT_DOWN)  //button released
        printf("BUTTON UP\n");
    else
        printf("BUTTON DOWN\n");  //button clicked
    
    //redraw the scene after mouse click
    glutPostRedisplay();
}


void read_data(string* data, int &counter, Polygon* arrPoly, int i)
{
	int numvertices;
	istringstream(data[counter++]) >> numvertices;
	arrPoly[i] = Polygon(numvertices); 

	for(int j = 0; j < numvertices; j++)
	{
		double x, y;
		istringstream(data[counter++]) >> x >> y;
		arrPoly[i].X[j] = x;
		arrPoly[i].Y[j] = y;
	}
}

void read_file(string* data)
{
	int count = 0;
	string line;
    ifstream myfile;

    myfile.open("test_scene"); //or test_scene 
    while (getline(myfile, line))
    	if(line.length() > 0) //only get the lines with info
    		data[count++] = line;

    myfile.close();
}

void check()
{
	GLenum err = glGetError();
	if(err != GL_NO_ERROR)
	{
		printf("GLERROR: There was an error %s\n",gluErrorString(err) );
		exit(1);
	}
}

void draw_pix(int x, int y){
    glBegin(GL_POINTS);     
	    glColor3f(.2,.2,1.0);  		//choose color 
	    glVertex3f(x+.5,y+.5,0);	//plot point
    glEnd();
}


//called repeatedly when glut isn't doing anything else
void idle()
{
    //redraw the scene over and over again
	glutPostRedisplay();	
}

/*initialize gl stufff*/
void init()
{
    //set clear color (Default background to white)
	glClearColor(1.0,1.0,1.0,1.0);
    //checks for OpenGL errors
	check();
}

/*Gets called when display size changes, including initial craetion of the display*/
void reshape(int width, int height)
{
	/*set up projection matrix to define the view port*/
    //update the ne window width and height
	win_width = width;
	win_height = height;
    
    //creates a rendering area across the window
	glViewport(0,0,width,height);
    // up an orthogonal projection matrix so that
    // the pixel space is mapped to the grid space
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	glOrtho(0,grid_width,0,grid_height,-10,10);
    
    //clear the modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    //set pixel size based on width, if the aspect ratio
    //changes this hack won't work as well
    pixel_size = width/(double)grid_width;
    
    //set pixel size relative to the grid cell size
    glPointSize(pixel_size);
    //check for opengl errors
	check();
}

int main(int argc, char** argv)
{

	//the number of pixels in the grid
    grid_width = 100;
    grid_height = 100;
    view_width = 100;
    view_wmin = 0;
    view_hmin = 0;
    view_height = 100;

//FIXME UNCOMMENT!!
    cout << "Specify grid dimensions, Enter width(enter 100 for default: ";
    cin >> grid_width;
    cout << "and height(enter 100 for default): ";
    cin >> grid_height;

    cout << "Specify viewport maximums.(Make it <= grid width and height)"<<endl<<"Width: ";
    cin >> view_width;
    cout << "Height: ";
    cin >> view_height;
    cout << "Specify viewport minimums (Make it >= default 0)" <<endl<<"Width: ";
    cin >> view_wmin;
    cout << "Height: ";
    cin >> view_hmin;                            

    cout << "How much do you want to scale/zoom in/out every time? Enter 0.2 if you want default: ";
    cin >> userscale;

    cout << "Translation vector x (enter 1 for default): ";
    cin >> translx;

    cout << "Translation vector y (enter 1 for default): ";
    cin >> transly;

    cout << "Enter the angle of rotation (enter 45 for default): ";
    cin >> angle;

    // the size of pixels sets the inital window height and width
    // don't make the pixels too large or the screen size will be larger than
    //your display size
    pixel_size = 5;

    /*Window information*/
    win_height = grid_height*pixel_size;
    win_width = grid_width*pixel_size;

    glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	/*initialize variables, allocate memory, create buffers, etc. */
    //create window of size (win_width x win_height
    glutInitWindowSize(win_width,win_height);
    //windown title is "glut demo"
	glutCreateWindow("glut demo");

	glutDisplayFunc(display); //then write a display func
	glutIdleFunc(idle);

	glutMouseFunc(mouse);     //mouse button events
	glutMotionFunc(motion);   //mouse movement events
	glutKeyboardFunc(key);    //Keyboard events

	glutReshapeFunc(reshape);	//includng this rly makes the shapes

	init();				//make bg white not black
	glutMainLoop();
}


void freeobjects(Polygon* poly, int** buffer, string* data, int n)
{
	for(int i = 0; i < n; i++)
	{
		delete[] poly[i].X;
		delete[] poly[i].Y;
	}

	delete [] poly;
	for(int i = 0; i < grid_width; i++)
		delete [] buffer[i];
	delete [] buffer;
	delete [] data;
}