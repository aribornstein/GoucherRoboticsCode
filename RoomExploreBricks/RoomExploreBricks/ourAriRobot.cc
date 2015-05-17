#include <iostream>
#include <fstream>
#include <libplayerc++/playerc++.h>
#include <GL/glut.h>
#include <algorithm>
#include <vector>
#include <math.h>

/* 
 * Some basic constants.  MAX_BALLS is rather meaningless at this point.
 * ESC is the ASCII value of the Esc key.  ELASTICITY is used to define
 * the elasticity of collisions.  It may range between 1.0 (completely
 * elastic) to 0.0 (completely inelastic).	VELOCITY_SCALE is used to scale
 * velocity to a reasonable value on fast machines.
 */

#define ESC 0x1b
#define VELOCITY_SCALE 0.33
#define LP_SIZE 654
#define LASER_Y_OFFSET_MM 45
#define ROBOT_RADIUS 165
#define WINDOW_SCALE 0.2
#define OBJECT_MAX 8	
#define BRICK_WIDTH 100 
#define BRICK_HEIGHT 50
#define BRICK_DIAG 120
#define FLOOR_TILE 30.5

using namespace std;
using namespace PlayerCc;

// TODO: PS
PlayerClient robot("anna.goucher.edu"); // rename to Anna PS
Position2dProxy pp(&robot,0);	
BumperProxy	bp(&robot,0);	
LaserProxy 	lp(&robot,0);	// only need for Anna PS

// TODO: fix here.
PlayerClient robot2("gumstix.goucher.edu");
Position2dProxy ottopp(&robot2,0);	
BumperProxy	ottobp(&robot2,0);	

//playerc_position2d_t *position; 

double lpReplacement[LP_SIZE];	// PS This array is our hack to fix robot.lp error

// TODO: Move me to a method or something
double lpBufferCopy[LP_SIZE];

/* 
 * Basic data structures for the simulation.  We would be better off doing
 * this in C++ and using proper classes.
 */

typedef struct Color
{
	GLdouble r, g, b;
} Color;

typedef pair <double, double> Pair;

inline bool less_than_second( const Pair& b1, const Pair& b2 ){
   return b1.second < b2.second;
}

/* This should really be extended to three dimensions */

typedef struct Vector2
{
	GLdouble x, y;
} Vector2;

/*
 * Most of these are self-explanatory.	handle is the display list necessary
 * for rendering a ball.
 */

typedef struct Ball //this is anna
{
	Vector2 position;	/* Ok, so it's not really a vector.	 Sue me. */ //I will //AB
	Vector2 velocity;
	GLdouble radius;
	GLdouble mass;
	Color color;
	GLuint handle;
} Ball;

typedef struct Room
{
	int height;
	int width;
} Room;

/***********************************************************************
 * Prototypes for the basic OpenGL functions.
 ***********************************************************************/

void display(void);
void init(void);
void reshape(int w, int h);
void idle(void);
void keyboard(unsigned char key, int x, int y);

Room getRoom(void);
void traverseBricks();

int range[66];

/* I created variables for anna and otto. Do we need that? Can we just use
   one copy and call them robotMove, etc etc. PS */
void annaMove(double amount);
void annaTurn(double amount);
double annaTheta;

double getNormalLength(double objX, double objY, double ottoX, double ottoY);
double getTheta(double objX, double objY, double ottoX, double ottoY);

//BlobfinderProxy bfp(&robot,0);

Room world;
Room annaPos;
//Room ottoPos; // PS
//Room ballsPos;

Ball anna;
//Ball otto; //PS
Ball brick[OBJECT_MAX];
Ball balls[OBJECT_MAX];

int pinkBricks[OBJECT_MAX];
int brickis = 0;
int objectCount = 0;

void placeBalls()
{
	anna.position.x = annaPos.width;
	anna.position.y = annaPos.height;
	anna.velocity.x = 0;
	anna.velocity.y = 0;
	//anna.velocity = scalarProduct(VELOCITY_SCALE, anna.velocity);
	//anna.position = scalarProduct(40.0, anna.

	// PS
	//otto.position.x = ottoPos.width;
	//otto.position.y = ottoPos.height;
	//otto.velocity.x = 0;
	//otto.velocity.y = 0;
}

void initBalls(void)
{
	GLUquadricObj *qobj;  /* Need this to generate the spheres (disks). */
	/* Compute starting positions and velocities for the balls. */

	placeBalls();

	anna.radius = ROBOT_RADIUS;
	//anna.mass = 2.0;
	anna.color.r = 1.0;
	anna.color.g = 0.0;
	anna.color.b = 0.0;
	/* Create display lists for each ball. */
	anna.handle = glGenLists(1);
	qobj = gluNewQuadric();
	glNewList(anna.handle, GL_COMPILE);
	gluDisk(qobj, 0.0, anna.radius, 72, 1);
	glColor3ub(0, 0, 255);
	gluPartialDisk(qobj, 0.0, anna.radius, 72, 1, 338, 45);
	glEndList();

	for (int i = 0; i < objectCount; i++) {
		balls[i].handle = glGenLists(1);
	
		// check if object is brick sized
		if (balls[i].radius > 20 && balls[i].radius < 80) {
			balls[i].color.r = 0.0;
			balls[i].color.g = 1.0;
			balls[i].color.b = 0.0;
			
			brick[brickis].position = balls[i].position;
			brick[brickis].position.x = world.width - brick[i].position.x;
			brick[brickis].position.x += BRICK_WIDTH / 2;
			brick[brickis].position.y += BRICK_HEIGHT / 2;
			
			brick[brickis].radius = balls[i].radius;	
			brickis++;			
		}
		// else its a robot.
		else {
			glNewList(balls[i].handle, GL_COMPILE);
			balls[i].color.r = 1.0;
			balls[i].color.g = 1.0;
			balls[i].color.b = 0.0;
			balls[i].radius = ROBOT_RADIUS; // Keep this?
			gluDisk(qobj, 0.0, balls[i].radius, 72, 1);
			glEndList();
		}
	}
	cout << "Brick count = " << brickis << endl;
}



/***********************************************************************
 * OpenGL functions.
 ***********************************************************************/

/***********************************************************************
 * Recall, this will do our rendering for us.  It is called following
 * each simulation step in order to update the window.
 ***********************************************************************/

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	/* Render balls.  why is this a mirror image*/
	glColor3f(anna.color.r, anna.color.g, anna.color.b);
	glPushMatrix();
	// manually flip x axis for robot. PS
	//glTranslatef(world.width - anna.position.x, anna.position.y, 0.0);
	glTranslatef(anna.position.x, anna.position.y, 0.0);
	glCallList(anna.handle);
	glPopMatrix();
	
	for (int i = 0; i < objectCount; i++) {
		/* Render balls.  why is this a mirror image*/
		glColor3f(balls[i].color.r, balls[i].color.g, balls[i].color.b);
		glPushMatrix();
		// manually flip x axis for robot. PS
		//glTranslatef(world.width - anna.position.x, anna.position.y, 0.0);
		glTranslatef(world.width - balls[i].position.x, balls[i].position.y, 0.0);
		glCallList(balls[i].handle);
		glPopMatrix();
	}

	//AB
	//glColor3f(otto.color.r, otto.color.g, otto.color.b);
	//glPushMatrix();
	// manually flip x axis for robot. PS
	//glTranslatef(world.width - otto.position.x, otto.position.y, 0.0);
	//glTranslatef(otto.position.x, otto.position.y, 0.0);
	//glCallList(otto.handle);
	//glPopMatrix();

	for(int i =0; i < brickis; i++){
		glBegin(GL_QUADS);
	
		if(pinkBricks[i]){
			glColor3f(2.6, 1.0, 1.5);}
		else{
			glColor3f(1.0, 0.0, 0.0);}
	
		// TODO: Fix orientation of brick
		// if radius > < something then one way
		// else other way
		
		cout << "brickis = " << i << endl;
		cout << "\t  @ " << brick[i].position.x << " " << brick[i].position.y << endl;
		glBegin(GL_QUADS);
		glVertex2f(brick[i].position.x, brick[i].position.y); // bottom left
		glVertex2f(brick[i].position.x, brick[i].position.y + BRICK_HEIGHT); // top left
		glVertex2f(brick[i].position.x + BRICK_WIDTH, brick[i].position.y + BRICK_HEIGHT); // top right
		glVertex2f(brick[i].position.x + BRICK_WIDTH, brick[i].position.y); // bottom left
		glEnd();
	}

	/* Swap buffers, for smooth animation.	This will also flush the
	 * pipeline.
	 */

	glutSwapBuffers();
}

void init(void)
{
	glClearColor (0.0, 0.0, 0.0, 0.0);
	glShadeModel (GL_FLAT);	  /* Probably unnecessary. */

	initBalls();
	traverseBricks();
}

/***********************************************************************
 * Hitting the Esc key will exit the program.
 ***********************************************************************/

// This doesn't work. PS
void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
		case ESC:
			break;
	}
}

/***********************************************************************
 * Set the basic world coordinates to screen coordinates mapping.
 ***********************************************************************/

void reshape(int w, int h)
{
	/* Probably needs to be fixed. */

	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, world.width, 0, world.height, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void idle(void)
{
	//anna.position.y += anna.velocity.y;
	//anna.position.x += anna.velocity.x;
	//anna.velocity.y = 0;
	//anna.velocity.x = 0;
	//glutPostRedisplay();
}

/************************************************************
 * anna functions
 *************************************************************/

void annaMove(double amount)
{
	// convert from mm to cm
	//amount /= 10;
	
	// If the amount is positive move forward, otherwise move backward
	int direction =1;
	if (amount<0)
	{
		direction =-1;
		amount = -amount;
	}

	//instance variables
	robot.Read();//first read
	double x = pp.GetXPos();// x position from first read
	double y = pp.GetYPos();// y position from first read
	double distance =0;
	double posDiff =0;

	//begin move process move
	pp.SetSpeed(direction*0.2,0.0);
	for(;;)
	{
		robot.Read();
		double x1 = pp.GetXPos();// x position from new read
		double y1 = pp.GetYPos();// y position from new read
		posDiff=sqrt(pow((x-x1),2)+pow((y-y1),2));
		distance = distance+ posDiff;
		x=x1; y=y1; // prepare for next pass
		if (distance> amount/145)return;
	}

	// stop when moving distance reached
	pp.SetSpeed(0.0,0.0);
}
void ottoMove(double amount)
{
	// convert from mm to cm
	amount /= 10;
	
	// If the amount is positive move forward, otherwise move backward
	int direction =1;
	if (amount<0)
	{
		direction =-1;
		amount = -amount;
	}

	//instance variables
	robot2.Read();//first read
	double x = ottopp.GetXPos();// x position from first read
	double y = ottopp.GetYPos();// y position from first read
	double distance =0;
	double posDiff =0;

	//begin move process move
	ottopp.SetSpeed(direction*0.2,0.0);
	for(;;)
	{
		robot2.Read();
		double x1 = ottopp.GetXPos();// x position from new read
		double y1 = ottopp.GetYPos();// y position from new read
		posDiff=sqrt(pow((x-x1),2)+pow((y-y1),2));
		distance = distance+ posDiff;
		x=x1; y=y1; // prepare for next pass
		if (distance> amount/145)return;
	}

	// stop when moving distance reached
	ottopp.SetSpeed(0.0,0.0);
}
/*
 * Get laser readings in externally then reads data from file to lpReplacement.
 * This is a work around for empty LP readings.
 */
void laserRead()
{
	// Delete old readings
	system("rm lpreadings.csv");

	// Get new reading
	system("./laserRead ");

	// Read from file
	ifstream inFile ("lpreadings.csv");
	if (inFile.is_open())
	{
		string line;
		int i = 0;
		while (inFile.good() && i < LP_SIZE)
		{
			// makesure this is secure
			getline (inFile,line);
			lpReplacement[i] = atof(line.c_str());
			i++;
		}
		inFile.close();
	}
}

void annaTurn(double amount)
{
	//constants
	const double absloute_amt =fabs(amount);
	//instance variables
	robot.Read();//first read
	double yaw1= pp.GetYaw();//Yaw from first read
	double yawDiff=0; 
	double angle = 0;
	
	//set speed of robot
	if (amount >= 0) {
		pp.SetSpeed(0.0,-.2);
	}
	else {
		pp.SetSpeed(0.0,.2);
	}
	//Implement the turn
	for(;;)
	{
		robot.Read();//new read
		double yaw2 = pp.GetYaw();//new read
		
		if (yaw1*yaw2 <0)//if yaws not both negitave or both positive
		{
			yawDiff = fabs(yaw1)+fabs(yaw2);
			if (yawDiff>M_PI) {
				yawDiff= 2*M_PI-yawDiff;
			}
		}
		else yawDiff = fabs(yaw1 - yaw2);{
			angle = angle + yawDiff; //add compenstion
		}
		
		if (angle > absloute_amt/2.6){
			return;//check if the turn is made if so exit
		}
		else {
			yaw1=yaw2;//else prepare for next pass
		}
	}
	//turn has been completed stop the robit
	pp.SetSpeed(0.0,0.0);
}

void ottoTurn(double amount)
{
	//constants
	const double absloute_amt =fabs(amount);
	//instance variables
	robot2.Read();//first read
	double yaw1= ottopp.GetYaw();//Yaw from first read
	double yawDiff=0; 
	double angle = 0;
	
	//set speed of robot
	if (amount >= 0) {
		ottopp.SetSpeed(0.0,-.2);
	}
	else {
		ottopp.SetSpeed(0.0,.2);
	}
	//Implement the turn
	for(;;)
	{
		robot2.Read();//new read
		double yaw2 = ottopp.GetYaw();//new read
		
		if (yaw1*yaw2 <0)//if yaws not both negitave or both positive
		{
			yawDiff = fabs(yaw1)+fabs(yaw2);
			if (yawDiff>M_PI) {
				yawDiff= 2*M_PI-yawDiff;
			}
		}
		else yawDiff = fabs(yaw1 - yaw2);{
			angle = angle + yawDiff; //add compenstion
		}
		
		if (angle > absloute_amt/2.6){
			return;//check if the turn is made if so exit
		}
		else {
			yaw1=yaw2;//else prepare for next pass
		}
	}
	//turn has been completed stop the robit
	ottopp.SetSpeed(0.0,0.0);
}

/************************************************************
 * Generic Robot functions
 *************************************************************/

void robotMove(double amount, Position2dProxy proxy, PlayerClient robit)
{
	// If the amount is positive move forward, otherwise move backward
	int direction =1;
	if (amount<0)
	{
		direction =-1;
		amount = -amount;
	}

	//instance variables
	robit.Read();//first read
	double x = proxy.GetXPos();// x position from first read
	double y = proxy.GetYPos();// y position from first read
	double distance =0;
	double posDiff =0;

	//begin move process move
	proxy.SetSpeed(direction*0.2,0.0);
	for(;;)
	{
		robit.Read();
		double x1 = proxy.GetXPos();// x position from new read
		double y1 = proxy.GetYPos();// y position from new read
		posDiff=sqrt(pow((x-x1),2)+pow((y-y1),2));
		distance = distance+ posDiff;
		x=x1; y=y1; // prepare for next pass
		if (distance> amount/145)return;
	}

	// stop when moving distance reached
	proxy.SetSpeed(0.0,0.0);
}

void robotTurn(double amount,Position2dProxy proxy, PlayerClient robit)
{
	//constants
	const double absloute_amt =fabs(amount);
	//instance variables
	robit.Read();//first read
	double yaw1= proxy.GetYaw();//Yaw from first read
	double yawDiff=0; 
	double angle = 0;
	
	//set speed of robot
	if (amount >= 0) {
		proxy.SetSpeed(0.0,-.2);
	}
	else {
		proxy.SetSpeed(0.0,.2);
	}
	//Implement the turn
	for(;;)
	{
		robit.Read();//new read
		double yaw2 = proxy.GetYaw();//new read
		
		if (yaw1*yaw2 <0)//if yaws not both negitave or both positive
		{
			yawDiff = fabs(yaw1)+fabs(yaw2);
			if (yawDiff>M_PI) {
				yawDiff= 2*M_PI-yawDiff;
			}
		}
		else yawDiff = fabs(yaw1 - yaw2);{
			angle = angle + yawDiff; //add compenstion
		}
		
		if (angle > absloute_amt/2.6){
			return;//check if the turn is made if so exit
		}
		else {
			yaw1=yaw2;//else prepare for next pass
		}
	}
	//turn has been completed stop the robit
	proxy.SetSpeed(0.0,0.0);
}

/*
 * Return an array that have the readings at 0 deg, 90 deg,
 * and 180 deg from the horizon. Values are an average of 3 neighbors.
 * Units are in mm.
 */
 // TODO: Fix me, signature method.
double * getFrontandSides(double array[])
{
	laserRead(); 

	int leftIndex =(LP_SIZE *5)/6;
	double leftValue = ((lpReplacement[leftIndex-1]
				+ lpReplacement[leftIndex]
				+ lpReplacement[leftIndex+1])/3)* 1000;

	int rightIndex = LP_SIZE/6;
	double rightValue =((lpReplacement[rightIndex-1]
				+ lpReplacement[rightIndex]
				+ lpReplacement[rightIndex+1])/3)* 1000;

	int frontIndex = LP_SIZE/2;
	double frontValue =((lpReplacement[frontIndex-1]
				+ lpReplacement[frontIndex]
				+ lpReplacement[frontIndex+1])/3)* 1000;

	array[0] = frontValue - LASER_Y_OFFSET_MM;
	array[1] = rightValue;
	// array[2] = behindValue. Not measurable in this method.
	array[3] = leftValue;

	return array;
}

/*
 * Returns the distance of the wall directly infront of the robot.
 * Units are mm.
 */
double getFront()
{
	laserRead();

	int frontIndex = LP_SIZE/2;
	double frontValue =((lpReplacement[frontIndex-1]
				+ lpReplacement[frontIndex]
				+ lpReplacement[frontIndex+1])/3)* 1000;
	return frontValue - LASER_Y_OFFSET_MM;
}

/* Returns the stdDev of two numbers. */
double sampleStdDev(double a, double b)
{
	double mean = (a + b) / 2;
	double stddev = sqrt( pow(a - mean, 2) + pow(b - mean, 2) );
	return stddev;
}

/*
 * WORK IN PROGRESS
 *
 * Find spikes in the stdDev of two data sets and returns the middle.
 * attempting to find a robot's center as an index of LP.
 */
 //int list findSpikes
void findSpikes(double* data, int* spikesList, int *spikeIndex)
{
	// TODO: Error handling when < or > two spikes are detected.

	for (int i = 0; i < LP_SIZE -1; i++)
	{
		if (data[i] < 4.1 && data[i+1] < 4.1){
			double diff =  data[i] - data[i+1];
			if (abs(diff) >= .1 && abs(diff) < 2)
			{	
				spikesList[*spikeIndex] = i;
				*spikeIndex = *spikeIndex + 1;
				cout << "found spike @" << i << endl;
			}
		}
		else cout <<"blip"<<endl;
	}
}

/* Finds the width of an object give the bounding indexes and distance. */
double findObjectWidth(int leftIndex, int rightIndex)
{
	/* Note: I move the leftIndex and rightIndex by one. 
	 * This is to avoid the error of scans over lapping and giving us weird
	 * distances for the edges. PROPER FIX NEEDED */
	 
	//double avgDistance = 0;
	//for (int i = leftIndex; i <= rightIndex; i++) {
	//	avgDistance += lpBufferCopy[i];
	//}
	//avgDistance /= (leftIndex - rightIndex);
	//double leftDistance = rightDistance = avgDistance;
	 
	double leftDistance = lpBufferCopy[leftIndex + 1] * 1000; // PS
	double rightDistance = lpBufferCopy[rightIndex - 1] * 1000; // PS
	
	double theta = abs(leftIndex - rightIndex) / (654.0 / 270.0);
	theta = theta * (M_PI / 180); // deg to rad

	double objectWidth = pow(leftDistance, 2) + pow(rightDistance, 2);
	objectWidth -= (2 * leftDistance * rightDistance * cos(theta));
	objectWidth = sqrt(objectWidth);
	// objectWidth; // TODO: Fix units, am I required?

	
	//cout << "Width Debug values" << endl;
	//cout << "\tLeft = " << leftDistance << " @ " << leftIndex << endl;
	//cout << "\tRight = " << rightDistance << " @ " << rightIndex << endl;
	//for (int i = leftIndex; i <= rightIndex; i++) {
	//	cout << "\tValue at " << i << " = " << lpBufferCopy[i] << endl; 
	//}
	//cout << "\tTheta = " << theta << endl;
	//cout << "\tWidth = " << objectWidth << endl;
	//*/
	
	return objectWidth;
}

/* 
 * WORK IN PROGRESS
 *
 * Rename? Perhaps findRobotCenter or findOttoCenter would be better name.
 * Gets two readings, the first with an empty room, the second with an obj in it.
 * Runs stdDev on all readings and then passes data to findSpike method.
 */
 // calcObjectsCenters
 

void scanRoom(int *spikeList, int * spikeIndex)
{
	// First read
	laserRead();

	double firstReading[LP_SIZE];	
	for (int i = 0; i < LP_SIZE; i++)
	{
		firstReading[i] = lpReplacement[i];
	}

	// Wait for key press
	cout << "First reading done. Insert robot then press any key." << endl;
	string input;
	getline(cin, input);
	cout << "Input received" << endl;

	// Second Read
	laserRead();

	double secondReading[LP_SIZE];
	for (int i = 0; i < LP_SIZE; i++)
	{
		secondReading[i] = lpReplacement[i];
		// save this buffer, it holds the distance to Otto
		lpBufferCopy[i] = lpReplacement[i];
	}

	// Calculate stdDev between readings, save to array, and write to file.
	ofstream outFile ("stddev.csv") ;

	double stdDevArray[LP_SIZE];
	for (int i = 0; i < LP_SIZE; i++)
	{
		stdDevArray[i] = sampleStdDev(firstReading[i], secondReading[i]);
		outFile << firstReading[i] << "\t" << secondReading[i] 
			<< "\t" << stdDevArray[i] << endl;
	}

	outFile.close();
	//check this
	findSpikes(stdDevArray, spikeList, spikeIndex);
}

void findObjectCenters(int* spikesList, int * result, int *spikeIndex)
{	
	int index= *spikeIndex;
	
	if (index % 2 == 0 && index > 0){
		// success
		// int objectCenters[(index) / 2]; 
		
		for (int i = 0; i < index; i+=2) {
			result[i/2] = (spikesList[i] + spikesList[i+1]) / 2;
		}
	}
	else {
		cout << "Error finding centers" << endl;
		pp.SetSpeed(0.0,0.0);
		exit (-1);
		int error[] = {-1};
		result = error; // fix me later
	}
	*spikeIndex = index;
}

// calcObjectPos
Room findObjectPosition(Room origin, int objectCenter)
{
	Room object;

	// Couldn't find Otto. Place him ontop of Anna for now.
	if (objectCenter == -1)
	{
		object.width = origin.width;
		object.height = origin.height;
	}
	else
	{
		// Center == the index of the center of LP
		double center = LP_SIZE / 2.0;

		// Find angle between object and the front of the origin.
		double theta = (((double) objectCenter - center) / (double) LP_SIZE) * 270;
		theta *= (M_PI / 180);

		// Converting from meters to millimeters.
		double distanceToobject = lpBufferCopy[objectCenter] * 1000;

		object.width = origin.width + (int) (distanceToobject * sin(theta));
		object.height = origin.height + (int) (distanceToobject * cos(theta));
	}
	
	return object;
}

/*
 * Calculates the size of the room and the robots positon within it.
 * Units are in mm?
 */
Room getRoom()
{
	// Get front, left, and right walls.
	double wall[4];
	
	float annaDistance = FLOOR_TILE * 3;
	annaMove(annaDistance);
	pp.SetSpeed(0.0,0.0);
	//wait one second
	sleep(1);
	getFrontandSides(wall);
	
	int spikesList [(OBJECT_MAX + 1)*2];
	int spikeIndex = 0; 
	cout << "about to call scanRoom" << endl;
	sleep(1);
	scanRoom(spikesList, &spikeIndex);
	
    int objectCenters[spikeIndex / 2];
	findObjectCenters(spikesList, objectCenters, &spikeIndex);
	objectCount = spikeIndex / 2;

	// Turn 180 deg
	annaTurn(M_PI);
	annaTurn(M_PI/32);
	pp.SetSpeed(0.0,0.0);
    sleep(1);
	// Get back wall.
	wall[2] = getFront();

	// Complete circle by turning 180 deg again
	annaMove(annaDistance * .9);
	annaTurn(M_PI);
	annaTurn(M_PI/32); // padding
	pp.SetSpeed(0.0,0.0);

	// The H,W of the room the robots are in.
	Room robotRoom;
	robotRoom.height = (int) (wall[0] + wall[2]);
	robotRoom.width = (int) (wall[1] + wall[3]);
	cout << "Room Width:" << robotRoom.width << " mm, Room Height:" 
		 << robotRoom.height << " mm." << endl;
	cout << "Room Area:" << robotRoom.height * robotRoom.width 
		 << " mm^2." << endl;

	// The X,Y of Anna within her room
	Room annaPosition; // get rid of annaPosition and just use global annaPos?
	annaPosition.height = (int) (robotRoom.height - wall[0]);
	annaPosition.width = (int) (robotRoom.width - wall[1]);
	cout << "Anna X:"<< annaPosition.width << " mm, Anna Y:"
		 << annaPosition.height << " mm." << endl;

	for (int i = 0; i < objectCount; i++){
		Room objectRoom = findObjectPosition(annaPosition, objectCenters[i]);
		balls[i].position.x = objectRoom.width + 650; //600 hackish offset so bricks dont go out of bounds
		balls[i].position.y = objectRoom.height;
		balls[i].radius = findObjectWidth(spikesList[i*2], spikesList[i*2 + 1]);
		balls[i].radius /= 2; // convert diameter to radius
		cout << "Object " << i << " radius: " << balls[i].radius << endl;
	}
	
	//cout << balls[0].position.x << " " << balls[0].position.y << endl;	
	
	//set this to the global variable
	annaPos = annaPosition;
	//ottoPos = ottoPosition;
	
	return robotRoom;
}

double getNormalLength(double objX, double objY, double ottoX, double ottoY){
	double normalLength = sqrt((objX*objX) + (objY*objY)) * sqrt((ottoX*ottoX) + (ottoY*ottoY));
	return normalLength;
}

double getTheta(double objX, double objY, double ottoX, double ottoY){
	double dotProduct = (objX * ottoX) + (objY * ottoY);
	double normalLength = sqrt((objX*objX) + (objY*objY)) * sqrt((ottoX*ottoX) + (ottoY*ottoY));
	double turn = acos(dotProduct / normalLength);
	return turn;
}

void traverseBricks(){
if (brickis >= 2){	
		vector<Pair> pairs(brickis + 1);
		
		pairs[0] = Pair(world.width - 225, 225);
		
		for (int i = 0; i < brickis; i++) {
			pairs[i + 1] = Pair(brick[i].position.x, brick[i].position.y);
		}
		
		cout << "Pair Size: " << brickis + 1 << endl;
		
		sort (pairs.begin(), pairs.end(), less_than_second);
		
		cout << "Bricks" << endl;
		for (int i = 0; i < brickis + 1; i++){
			cout << pairs[i].first << " " << pairs[i].second << endl;
		}
		cout << endl;
		
		for (int i = 0; i < brickis; i++) {
			// Otto
			double a = 0;
			double b = 1;
			
			// Brick
			double c = pairs[i + 1].first - pairs[i].first;
			double d = pairs[i + 1].second - pairs[i].second;
			
			cout << a << " " << b << " " << c << " " << d << endl;
			
			double gNL = getNormalLength(a,b,c,d);
			cout << "gNL " << gNL << endl;
			
			double gTurn = getTheta(a,b,c,d);

			if (c < 0){ 
				gTurn *= -1;
			}
				
			cout << "Turn rad " << gTurn << endl;
			cout << "Turn deg " << gTurn * (180 / M_PI) << endl << endl;

			//to our future selves we reversed the polarity
			ottoTurn(gTurn);
			ottopp.SetSpeed(0.0,0.0);

			ottoMove(gNL);
			ottopp.SetSpeed(0.0, 0.0);
			
			ottoTurn(-gTurn);
			ottopp.SetSpeed(0.0, 0.0);
		}
	}
	else if (brickis > 0 && brickis < 2){
		//(a,b) = robot vector (c d) = brick vector
		double a = 0;
		double b = 1;
		// Otto positioned 5 cm from each wall.
		double c = balls[0].position.x - 225;
		double d = balls[0].position.y - 225;
		
		double dotProduct = (a*c) + (b*d);
		cout << "DP " << dotProduct << endl;
		 
		double normalLength = sqrt((c*c) + (d*d));
		cout << "NL " << normalLength << endl;
		
		double turn = acos(dotProduct / normalLength);
		cout << "Turn rad " << turn << endl;
		cout << "Turn deg " << turn * 180 / M_PI << endl;
		
		ottoTurn(-turn);
		ottopp.SetSpeed(0.0,0.0);
		/*
		ottoMove(normalLength + ROBOT_RADIUS);
		ottopp.SetSpeed(0.0, 0.0);
		*/
	}	
}

/* Here be dragons */
int main(int argc, char** argv)
{	 
	world = getRoom();
	// connect to otto here if he exists
	
	srand((unsigned int) time(NULL));
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	int windowWidth = (int) (world.width * WINDOW_SCALE);
	int windowHeight = (int) (world.height * WINDOW_SCALE);
	glutInitWindowSize (windowWidth, windowHeight);
	glutInitWindowPosition (100, 100);
	glutCreateWindow ("Anna's World");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);
	glutMainLoop();

	return 0;
}
