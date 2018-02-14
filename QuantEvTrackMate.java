import java.io.File;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;

import jsc.datastructures.*;
import jsc.tests.*;
import jsc.onesample.*;

public class QuantEv{

	public class CubicInterpolator
	{
		public double getValue (double[] p, double x) {
			return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
		}
	}

	// morphological mathematics
	// structuring element
	int[][] makeStrElt(int radiusXY,int radiusZ){
		int StrEltDimXY = 2*radiusXY+1,StrEltDimZ=2*radiusZ+1;
		int[][] StrElt = new int[StrEltDimZ][StrEltDimXY*StrEltDimXY];
		for(int dz=-radiusZ;dz<=radiusZ;dz++){
			for(int dy=-radiusXY;dy<=radiusXY;dy++){
				for(int dx=-radiusXY;dx<=radiusXY;dx++){
					if(Math.sqrt(Math.pow(dx,2d)+Math.pow(dy,2d)+Math.pow(dz,2d))<radiusXY){
						StrElt[dz+radiusZ][(dy+radiusXY)*StrEltDimXY+dx+radiusXY] = 1;
					}
				}
			}
		}
		return StrElt;
	}

	// dilation
	int[][] dilate(int[][] data,int[][] StrElt,int radiusXY,int radiusZ,int dimX,int dimY,int dimZ){
		int[][] result = new int[dimZ][dimX*dimY];
		int StrEltDimXY = 2*radiusXY+1,StrEltDimZ = 2*radiusZ+1;
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				for(int x=0;x<dimX;x++){
					if(data[z][y*dimX+x]>0){
						for(int w=0;w<StrEltDimZ;w++){
							for(int v=0;v<StrEltDimXY;v++){
								for(int u=0;u<StrEltDimXY;u++){
									if(((x+u-radiusXY)>=0)&&((x+u-radiusXY)<dimX)&&((y+v-radiusXY)>=0)&&((y+v-radiusXY)<dimY)&&((z+w-radiusZ)>=0)&&((z+w-radiusZ)<dimZ)){
										if(StrElt[w][v*StrEltDimXY+u]>0.0001){
											result[z+w-radiusZ][(y+v-radiusXY)*dimX+x+u-radiusXY] = 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		return result;
	}

	// compute cell border
	int[][] computePseudo3DCellBorder(int[][] mask,int dimX,int dimY,int dimZ){
		int strEltXY = 2,strEltZ = 1;
		if(dimZ>1){strEltZ=2;}
		int[][] strElt = makeStrElt(strEltXY, strEltZ);
		int[][] output = dilate(mask,strElt,strEltXY,strEltZ,dimX,dimY,dimZ);
		for(int z=0;z<dimZ;z++){
			for(int i=0;i<(dimX*dimY);i++){
				output[z][i] = output[z][i] - mask[z][i];
			}
		}

		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				if(mask[z][y*dimX]==1){output[z][y*dimX] = 1;}
				if(mask[z][y*dimX+dimX-1]==1){output[z][y*dimX+dimX-1] = 1;}
			}
			for(int x=0;x<dimX;x++){
				if(mask[z][x]==1){output[z][x] = 1;}
				if(mask[z][(dimY-1)*dimX+x]==1){output[z][(dimY-1)*dimX+x] = 1;}
			}
		}
		return output;
	}

	int[] extractXcellBorderCoordinates(int[][] cellBorder,int dimX,int dimY, int dimZ){
		int nbTrajectories=0;
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				for(int x=0;x<dimX;x++){
					if(cellBorder[z][y*dimX+x]>0){
						nbTrajectories++;
					}
				}
			}
		}
		int[] output = new int[nbTrajectories];
		int i=0;
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				for(int x=0;x<dimX;x++){
					if(cellBorder[z][y*dimX+x]>0){
						output[i] = x;
						i++;
					}
				}
			}
		}
		return output;
	}
	int[] extractYcellBorderCoordinates(int[][] cellBorder,int dimX,int dimY, int dimZ){
		int nbTrajectories=0;
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				for(int x=0;x<dimX;x++){
					if(cellBorder[z][y*dimX+x]>0){
						nbTrajectories++;
					}
				}
			}
		}
		int[] output = new int[nbTrajectories];
		int i=0;
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				for(int x=0;x<dimX;x++){
					if(cellBorder[z][y*dimX+x]>0){
						output[i] = y;
						i++;
					}
				}
			}
		}
		return output;
	}

	// compute cell border
	int[][] compute3DCellBorder(int[][] mask,int dimX,int dimY,int dimZ){
		int strEltXY = 2,strEltZ = 1;
		if(dimZ>1){strEltZ=2;}
		int[][] strElt = makeStrElt(strEltXY, strEltZ);
		int[][] output = dilate(mask,strElt,strEltXY,strEltZ,dimX,dimY,dimZ);
		for(int z=0;z<dimZ;z++){
			for(int i=0;i<(dimX*dimY);i++){
				output[z][i] = output[z][i] - mask[z][i];
			}
		}

		for(int i=0;i<(dimX*dimY);i++){
			if(mask[0][i]==1){output[0][i]=1;}
			if(mask[dimZ-1][i]==1){output[dimZ-1][i]=1;}
		}
		for(int z=1;z<(dimZ-1);z++){
			for(int y=0;y<dimY;y++){
				if(mask[z][y*dimX]==1){output[z][y*dimX] = 1;}
				if(mask[z][y*dimX+dimX-1]==1){output[z][y*dimX+dimX-1] = 1;}
			}
			for(int x=0;x<dimX;x++){
				if(mask[z][x]==1){output[z][x] = 1;}
				if(mask[z][(dimY-1)*dimX+x]==1){output[z][(dimY-1)*dimX+x] = 1;}
			}
		}

		return output;
	}

	int[][] computePseudo3DInnerBorder(int[][] mask,int dimX,int dimY,int dimZ){
		int strEltXY = 2,strEltZ = 1;
		if(dimZ>1){strEltZ=2;}
		int[][] strElt = makeStrElt(strEltXY, strEltZ);
		int[][] output = dilate(mask,strElt,strEltXY,strEltZ,dimX,dimY,dimZ);
		for(int z=0;z<dimZ;z++){
			for(int i=0;i<(dimX*dimY);i++){
				output[z][i] = output[z][i] - mask[z][i];
			}
		}
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				if(mask[z][y*dimX]==1){output[z][y*dimX] = 1;}
				if(mask[z][y*dimX+dimX-1]==1){output[z][y*dimX+dimX-1] = 1;}
			}
			for(int x=0;x<dimX;x++){
				if(mask[z][x]==1){output[z][x] = 1;}
				if(mask[z][(dimY-1)*dimX+x]==1){output[z][(dimY-1)*dimX+x] = 1;}
			}
		}
		
		return output;
	}
	
	int[][] compute3DInnerBorder(int[][] mask,int dimX,int dimY,int dimZ){
		int strEltXY = 2,strEltZ = 1;
		if(dimZ>1){strEltZ=2;}
		int[][] strElt = makeStrElt(strEltXY, strEltZ);
		int[][] output = dilate(mask,strElt,strEltXY,strEltZ,dimX,dimY,dimZ);
		for(int z=0;z<dimZ;z++){
			for(int i=0;i<(dimX*dimY);i++){
				output[z][i] = output[z][i] - mask[z][i];
			}
		}
		for(int i=0;i<(dimX*dimY);i++){
			if(mask[0][i]==1){output[0][i]=1;}
			if(mask[dimZ-1][i]==1){output[dimZ-1][i]=1;}
		}
		for(int z=1;z<(dimZ-1);z++){
			for(int y=0;y<dimY;y++){
				if(mask[z][y*dimX]==1){output[z][y*dimX] = 1;}
				if(mask[z][y*dimX+dimX-1]==1){output[z][y*dimX+dimX-1] = 1;}
			}
			for(int x=0;x<dimX;x++){
				if(mask[z][x]==1){output[z][x] = 1;}
				if(mask[z][(dimY-1)*dimX+x]==1){output[z][(dimY-1)*dimX+x] = 1;}
			}
		}
		
		return output;
	}

	// compute distance to cell border
	float computeDistanceToCellBorder(int referenceCenter1,int referenceCenter2,int x,int y,int z,int[][] cellBorder,int dimX,int dimY,int dimZ){
		boolean verticalLine=false;
		float a=0,b=0;
		if(x!=referenceCenter1){
			a=((float)(y-referenceCenter2)/(float)(referenceCenter1-x));
			b=(-referenceCenter2-a*referenceCenter1);
		}
		else{
			verticalLine=true;
		}

		float minOutsideDistance=10000;
		int outsideIdX=0,outsideIdY=0;

		for(int j=0;j<dimY;j++){
			for(int i=0;i<dimX;i++){
				if(cellBorder[z][j*dimX+i]>0){
					if((Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.))+Math.sqrt(Math.pow(x-referenceCenter1,2.)+Math.pow(y-referenceCenter2,2.)))<(Math.sqrt(Math.pow(i-referenceCenter1,2.)+Math.pow(j-referenceCenter2,2.)))+1.){
						float distance=0;
						if(verticalLine){
							if(i==x){
								distance=(float)(Math.sqrt(Math.pow((float)(j-referenceCenter2),2.)));
								if(distance<minOutsideDistance){
									minOutsideDistance = distance;
									outsideIdX = i;
									outsideIdY = j;
								}
							}
						}
						else{
							distance=(float)(Math.abs(a*i+j+b)/Math.sqrt(1+a*a));
							if(distance<minOutsideDistance){
								minOutsideDistance = distance;
								outsideIdX = i;
								outsideIdY = j;
							}
						}
					}
				}
			}
		}
		return (float)(Math.sqrt(Math.pow(referenceCenter1-outsideIdX,2.)+Math.pow(referenceCenter2-outsideIdY,2.)));
	}

	
	float computeDistanceToCellBorder(int referenceCenter1,int referenceCenter2,int x,int y,int z,int[][] cellBorder,int[][] innerBorder,int dimX,int dimY,int dimZ,int[][] insideCoordX,int[][] insideCoordY){
		boolean verticalLine=false;
		float a=0,b=0;
		if(x!=referenceCenter1){
			a=((float)(y-referenceCenter2)/(float)(referenceCenter1-x));
			b=(-referenceCenter2-a*referenceCenter1);
		}
		else{
			verticalLine=true;
		}

		float minOutsideDistance=10000,minInsideDistance=10000;
		int outsideIdX=0,outsideIdY=0,insideIdX=0,insideIdY=0;

		for(int j=0;j<dimY;j++){
			for(int i=0;i<dimX;i++){
				if(innerBorder[z][j*dimX+i]>0){
					if((Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.))+Math.sqrt(Math.pow(i-referenceCenter1,2.)+Math.pow(j-referenceCenter2,2.)))<(Math.sqrt(Math.pow(x-referenceCenter1,2.)+Math.pow(y-referenceCenter2,2.)))+1.){
						float distance=0;
						if(verticalLine){
							if(i==x){
								distance=(float)(Math.sqrt(Math.pow((float)(j-referenceCenter2),2.)));
							}
						}
						else{
							distance=(float)(Math.abs(a*i+j+b)/Math.sqrt(1+a*a));
						}
						if(distance<minInsideDistance){
							minInsideDistance = distance;
							insideIdX = i;
							insideIdY = j;
						}
					}
				}
				if(cellBorder[z][j*dimX+i]>0){
					if((Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.))+Math.sqrt(Math.pow(x-referenceCenter1,2.)+Math.pow(y-referenceCenter2,2.)))<(Math.sqrt(Math.pow(i-referenceCenter1,2.)+Math.pow(j-referenceCenter2,2.)))+1.){
						float distance=0;
						if(verticalLine){
							if(i==x){
								distance=(float)(Math.sqrt(Math.pow((float)(j-referenceCenter2),2.)));
							}
						}
						else{
							distance=(float)(Math.abs(a*i+j+b)/Math.sqrt(1+a*a));
						}
						if(distance<minOutsideDistance){
							minOutsideDistance = distance;
							outsideIdX = i;
							outsideIdY = j;
						}

					}
				}
			}
		}
		insideCoordX[z][y*dimX+x] = insideIdX;
		insideCoordY[z][y*dimX+x] = insideIdY;
		return (float)(Math.sqrt(Math.pow(insideIdX-outsideIdX,2.)+Math.pow(insideIdY-outsideIdY,2.)));
	}

	float[] computeDistanceToCellBorder(int[] cellBorderX,int[] cellBorderY,int nbTrajectories,int x, int y){
		float[] output = new float[nbTrajectories*3];
		float angleUnit=(float)(2*Math.PI/(float)nbTrajectories);
		for(int u=0;u<cellBorderX.length;u++){
			float theta=(float)Math.atan2(y-cellBorderY[u],cellBorderX[u]-x);
			if(theta<0){theta+=(2*Math.PI);}
			if(theta>(2*Math.PI)){theta-=(2*Math.PI);}
			output[(int)(theta/angleUnit)] = theta;
			output[(int)(theta/angleUnit)+nbTrajectories] = cellBorderX[u];
			output[(int)(theta/angleUnit)+2*nbTrajectories] = cellBorderY[u];
		}
		for(int i=0;i<nbTrajectories;i++){
			if((output[i]<0.001)&&(output[i+nbTrajectories]<0.001)&&(output[i+2*nbTrajectories]<0.001)){
				int previousIndex=-1,nextIndex=-1;
				boolean previousOver=false,nextOver=false;
				for(int u=i;u>=0;u--){
					if((output[u+nbTrajectories]>0.001)&&(output[u+2*nbTrajectories]>0.001)&&(!previousOver)){
						previousIndex=u;
						previousOver=true;
					}
				}
				for(int u=i;u<nbTrajectories;u++){
					if((output[u+nbTrajectories]>0.001)&&(output[u+2*nbTrajectories]>0.001)&&(!nextOver)){
						nextIndex=u;
						nextOver=true;
					}
				}
				output[i] = angleUnit*(float)i;
				if((previousIndex>-1)&&(nextIndex>-1)){
					output[i+nbTrajectories] = (output[previousIndex+nbTrajectories]+output[nextIndex+nbTrajectories])/2;
					output[i+2*nbTrajectories] = (output[previousIndex+2*nbTrajectories]+output[nextIndex+2*nbTrajectories])/2;
				}
				else{
					if(previousIndex>-1){
						output[i+nbTrajectories] = output[previousIndex+nbTrajectories];
						output[i+2*nbTrajectories] = output[previousIndex+2*nbTrajectories];
					}
					else{
						if(nextIndex>-1){
							output[i+nbTrajectories] = output[nextIndex+nbTrajectories];
							output[i+2*nbTrajectories] = output[nextIndex+2*nbTrajectories];
						}
					}
				}
			}
		}
		return output;
	}

	float computeSphericalDistanceToCellBorder(int referenceCenter1,int referenceCenter2,int referenceCenter3,int x,int y,int z,int[][] cellBorder,int dimX,int dimY,int dimZ){

		float minOutsideDistance=10000;
		int outsideIdX=0,outsideIdY=0,outsideIdZ=0;

		for(int k=0;k<dimZ;k++){
			for(int j=0;j<dimY;j++){
				for(int i=0;i<dimX;i++){
					if(cellBorder[k][j*dimX+i]>0){
						float distance = (float)Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.)+Math.pow(z-k,2.));
						if(distance<minOutsideDistance){
							minOutsideDistance = distance;
							outsideIdX = i;
							outsideIdY = j;
							outsideIdZ = k;
						}
					}
				}
			}
		}
		return (float)(Math.sqrt(Math.pow(referenceCenter1-outsideIdX,2.)+Math.pow(referenceCenter2-outsideIdY,2.)+Math.pow(referenceCenter3-outsideIdZ,2.)));
	}

	
	float computeSphericalDistanceToCellBorder(int referenceCenter1,int referenceCenter2,int referenceCenter3,int x,int y,int z,int[][] cellBorder,int[][] innerBorder,int dimX,int dimY,int dimZ){
		float minInsideDistance=10000,minOutsideDistance=10000;
		int insideIdX=0,insideIdY=0,insideIdZ=0,outsideIdX=0,outsideIdY=0,outsideIdZ=0;

		for(int k=0;k<dimZ;k++){
			for(int j=0;j<dimY;j++){
				for(int i=0;i<dimX;i++){
					if(cellBorder[k][j*dimX+i]>0){
						float distance = (float)Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.)+Math.pow(z-k,2.));
						if(distance<minOutsideDistance){
							minOutsideDistance = distance;
							outsideIdX = i;
							outsideIdY = j;
							outsideIdZ = k;
						}
					}
				}
			}
		}

		for(int k=0;k<dimZ;k++){
			for(int j=0;j<dimY;j++){
				for(int i=0;i<dimX;i++){
					if(innerBorder[k][j*dimX+i]>0){
						float distance = (float)Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.)+Math.pow(z-k,2.));
						if(distance<minInsideDistance){
							minInsideDistance = distance;
							insideIdX = i;
							insideIdY = j;
							insideIdZ = k;
						}
					}
				}
			}
		}
		
		return (float)(Math.sqrt(Math.pow(insideIdX-outsideIdX,2.)+Math.pow(insideIdY-outsideIdY,2.)+Math.pow(insideIdZ-outsideIdZ,2.)));
	}
	

	float computeActualDistanceToCellBorder(int x,int y,int z,int[][] cellBorder,int dimX,int dimY,int dimZ){

		float minOutsideDistance=10000;

		for(int k=0;k<dimZ;k++){
			for(int j=0;j<dimY;j++){
				for(int i=0;i<dimX;i++){
					if(cellBorder[k][j*dimX+i]>0){
						float distance = (float)Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.)+Math.pow(z-k,2.));
						if(distance<minOutsideDistance){
							minOutsideDistance = distance;
						}
					}
				}
			}
		}
		return minOutsideDistance;
	}

	
	float computeActualDistanceToCellBorder(int x,int y,int z,int[][] cellBorder,int[][] innerBorder,int dimX,int dimY,int dimZ){
		float minDistance=10000;

		for(int k=0;k<dimZ;k++){
			for(int j=0;j<dimY;j++){
				for(int i=0;i<dimX;i++){
					if(cellBorder[k][j*dimX+i]>0){
						float distance = (float)Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.)+Math.pow(z-k,2.));
						if(distance<minDistance){
							minDistance = distance;
						}
					}
				}
			}
		}

		for(int k=0;k<dimZ;k++){
			for(int j=0;j<dimY;j++){
				for(int i=0;i<dimX;i++){
					if(innerBorder[k][j*dimX+i]>0){
						float distance = (float)Math.sqrt(Math.pow(x-i,2.)+Math.pow(y-j,2.)+Math.pow(z-k,2.));
						if(distance<minDistance){
							minDistance = distance;
						}
					}
				}
			}
		}
		
		return minDistance;
	}

	// extract event coordinates from images
	void extractCylindricalCoordinates(float[][] input,int x,int y,int z,float rad,float orientation,float segment,boolean considerIntensity,
			float[] radius,float[] angle,float[] depth,float[] intensity,float[] distanceToBorder,
			int dimX,int dimY,int dimZ,int i){

			if(input[z][y*dimX+x]>0){
				angle[i] = orientation;
				if(considerIntensity){intensity[i] = input[z][y*dimX+x];}
				else{intensity[i] = 1;}
				depth[i] = z;
				radius[i] = rad;
				distanceToBorder[i] = segment;
			}
	}


	// compute angle and distance
	void computeAngleAndDistance(float[] angle,float[] distance,float[] intensitySum,int[][] mask,float[] distanceForAngle,int xRef,int yRef,double angleUnit,int dimX,int dimY,int nbTrajectories){
		int i=0;
		for(int y=0;y<dimY;y++){
			for(int x=0;x<dimX;x++){
				if((mask[0][y*dimX+x]>0)&&(intensitySum[y*dimX+x]>0)){
					double theta = Math.atan2(yRef-y,x-xRef);
					if(theta<0){theta+=(2*Math.PI);}
					if(theta>(2*Math.PI)){theta-=(2*Math.PI);}
					angle[i] = (float)theta;
					CubicInterpolator interpolator = new CubicInterpolator();
					double[] coordinates = new double[4];
					int lowerUnit=(int)Math.floor(theta/angleUnit),lowerUnitM=lowerUnit-1,upperUnit=lowerUnit+1,upperUnitP=lowerUnit+2;
					if(lowerUnitM<0){lowerUnitM=nbTrajectories-1;}
					if(upperUnit>(nbTrajectories-1)){upperUnit -= nbTrajectories;}
					if(upperUnitP>(nbTrajectories-1)){upperUnitP -= nbTrajectories;}
					coordinates[0] = distanceForAngle[nbTrajectories+lowerUnitM];
					coordinates[1] = distanceForAngle[nbTrajectories+lowerUnit];
					coordinates[2] = distanceForAngle[nbTrajectories+upperUnit];
					coordinates[3] = distanceForAngle[nbTrajectories+upperUnitP];
					double lowerMultiplier = theta/(float)angleUnit-lowerUnit,xCoord=interpolator.getValue(coordinates,lowerMultiplier);
					coordinates[0] = distanceForAngle[2*nbTrajectories+lowerUnitM];
					coordinates[1] = distanceForAngle[2*nbTrajectories+lowerUnit];
					coordinates[2] = distanceForAngle[2*nbTrajectories+upperUnit];
					coordinates[3] = distanceForAngle[2*nbTrajectories+upperUnitP];
					double yCoord=interpolator.getValue(coordinates,lowerMultiplier);
					double currentDistance=Math.sqrt(Math.pow(xRef-xCoord,2.)+Math.pow(yRef-yCoord,2.) );
					distance[i] = (float)currentDistance;
					i++;
				}
			}
		}
	}

	// compute angle and distance
	void computeAngleAndDistance(float[] angle,float[] distance,float[][] intensity,int[][] mask,float[] distanceForAngle,int xRef,int yRef,double angleUnit,int dimX,int dimY,int dimZ,int nbTrajectories){
		int i=0;
		for(int z=0;z<dimZ;z++){
			for(int y=0;y<dimY;y++){
				for(int x=0;x<dimX;x++){
					if((mask[z][y*dimX+x]>0)&&(intensity[z][y*dimX+x]>0)){
						float theta = (float)Math.atan2(yRef-y,x-xRef);
						if(theta<0){theta+=(2*Math.PI);}
						if(theta>(2*Math.PI)){theta-=(2*Math.PI);}
						angle[i] = theta;
						int lowerUnit=(int)Math.floor(theta/angleUnit);
						float lowerMultiplier = theta/(float)angleUnit-lowerUnit;							
						double currentDistance=Math.sqrt(Math.pow(xRef-(distanceForAngle[nbTrajectories+lowerUnit]*lowerMultiplier+distanceForAngle[nbTrajectories+lowerUnit+1]*(1-lowerMultiplier)),2.)
								+Math.pow(xRef-(distanceForAngle[2*nbTrajectories+lowerUnit]*lowerMultiplier+distanceForAngle[2*nbTrajectories+lowerUnit+1]*(1-lowerMultiplier)),2.) );
						distance[i] = (float)currentDistance;
						i++;
					}
				}
			}
		}
	}

	// extract event coordinates from images
	void extractSphericalCoordinates(float[][] input,int x,int y,int z,float rad,float orientation,float segment,boolean considerIntensity,
			float[] radius,float[] angle1,float[] angle2,float[] intensity,float[] distanceToBorder,
			int dimX,int dimY,int dimZ,int i){

			angle1[i] = orientation;
			if(considerIntensity){intensity[i] = input[z][y*dimX+x];}
			else{intensity[i] = 1;}
			angle2[i] = (float)Math.acos(z/Math.sqrt(rad*rad+z*z));
			radius[i] = (float)Math.sqrt(rad*rad+z*z);
			distanceToBorder[i] = segment;
	}
	
	// extract event coordinates from images
	void extractDistanceCoordinates(float[][] input,int x,int y,int z,float segment,boolean considerIntensity,
			float[] distance,float[] intensity,float[] unity,int dimX,int dimY,int dimZ,int i){

			distance[i] = segment;
			if(considerIntensity){intensity[i] = input[z][y*dimX+x];}
			else{intensity[i] = 1;}
			unity[i] = 1;
	}

		
	// extract event coordinates from images
	void extractCartesianCoordinates(float[][] input,int x,int y,int z,float rad,float orientation,float segment,boolean considerIntensity,
			float[] xTab,float[] yTab,float[] depth,float[] intensity,float[] distanceToBorder,
			int dimX,int dimY,int dimZ,int i){

				xTab[i] = rad*(float)Math.cos(orientation);
				yTab[i] = rad*(float)Math.sin(orientation);
				if(considerIntensity){intensity[i] = input[z][y*dimX+x];}
				else{intensity[i] = 1;}
				depth[i] = z;
				distanceToBorder[i] = segment;
	}


	// compute bandwidth according to the rule of thumb
	float compute_rule_of_thumb_bandwidth(float[] input) {

		float mu=0;

		// Compute mean and standard deviation
		float mean = 0, std = 0;
		for(int i=0;i<input.length;i++){
			mean += input[i];
		}
		mean /= (float)(input.length-1);
		for(int i=0;i<input.length;i++){
			std += (float)(Math.pow(((double)input[i]-mean),2.));
		}
		std = (float)Math.sqrt((double)std/(double)(input.length-1));
		mu = (float)(Math.pow((4.0*Math.pow(std,5.))/(3.0*(double)input.length), 1.0/5.0));

		return mu;
	}

	// compute histogram for radius
	float[] computeRadiusHistogram(float[] input,float[] weight,float[] length,int nbBins){

		float[] histogram = new float[nbBins];
		int[] nbEltsPerBin = new int[nbBins];
		float maxValue=0,lengthMean=0,lengthVar=0,nbEvtsPerExp=0;
		boolean takeLengthIntoAccount=false;
		for(int i=0;i<input.length;i++){
			lengthMean += length[i];
			nbEvtsPerExp += 1.;
		}
		lengthMean /= nbEvtsPerExp;
		for(int i=0;i<input.length;i++){
			lengthVar += Math.pow(length[i]-lengthMean,2.);
		}
		lengthVar /= nbEvtsPerExp;
		if(lengthVar>1.){takeLengthIntoAccount=true;}

		if(!takeLengthIntoAccount){
			for(int i=0;i<input.length;i++){
				if(input[i]>maxValue){maxValue=input[i];}
			}
		}
		for(int i=0;i<input.length;i++){
			int currentBin;
			if(takeLengthIntoAccount){
				if((!Float.isNaN(input[i]))&&(!Float.isNaN(length[i]))){
					currentBin = Math.round((float)(nbBins-1)*input[i]/length[i]);
					if(!Float.isNaN(currentBin)){
						if(currentBin==nbBins){currentBin=0;}
						if(!Float.isNaN(weight[i])&&!Float.isNaN(length[i])){
							histogram[currentBin] += weight[i];
							nbEltsPerBin[currentBin]++;
						}
					}
				}
			}
			else{
				if(!Float.isNaN(input[i])){
					currentBin = Math.round((float)(nbBins-1)*input[i]/maxValue);
					if(!Float.isNaN(currentBin)){
						if(currentBin==nbBins){currentBin=0;}
						if(!Float.isNaN(weight[i])&&!Float.isNaN(length[i])&&(weight[i]<1000000)){
							histogram[currentBin] += weight[i];
							nbEltsPerBin[currentBin]++;
						}
					}
				}
			}
		}

		for(int i=0;i<nbBins;i++){
			if(nbEltsPerBin[i]>0){
				histogram[i] /= (float)nbEltsPerBin[i];
			}
		}
		return histogram;
	}


	// compute histogram for radius
	float[] computeDepthHistogram(float[] input,float[] weight,float[] length,int nbBins,float maxValue,float[] interpolatedHistogram){

		int actualNbBins;
		if(nbBins>(maxValue+1.)){actualNbBins = (int)Math.ceil(maxValue+1);}
		else{actualNbBins = nbBins;}

		float[] histogram = new float[actualNbBins],
				outputHistogram = new float[nbBins];
		int[] nbEltsPerBin = new int[nbBins];
		for(int i=0;i<input.length;i++){
			if(!Float.isNaN(input[i])){
				int currentBin = Math.round((float)(actualNbBins-1)*input[i]/maxValue);
				if(currentBin==actualNbBins){currentBin=0;}
				if(!Float.isNaN(weight[i])&&!Float.isNaN(length[i])&&(weight[i]<1000000)){
					histogram[currentBin] += weight[i];
					nbEltsPerBin[currentBin]++;
				}
			}
		}
		double totalSum=0.;
		for(int i=0;i<actualNbBins;i++){
			if(nbEltsPerBin[i]>0){
				histogram[i] /= (float)nbEltsPerBin[i];
			}
		}

		if(actualNbBins!=nbBins){
			int translationBin = (int)((float)nbBins/(2*(float)actualNbBins));
			for(int i=0;i<translationBin;i++){
				interpolatedHistogram[i] = histogram[0];
			}
			for(int i=0;i<nbBins;i++){
				outputHistogram[i] = histogram[(int)((float)i*(float)actualNbBins/(float)nbBins)];
				int lowerBound=(int)Math.floor((float)i*(float)actualNbBins/(float)nbBins),
						higherBound=(int)Math.ceil((float)i*(float)actualNbBins/(float)nbBins);
				if(i<(nbBins-translationBin)){
					if(lowerBound==higherBound){
						interpolatedHistogram[i+translationBin] = histogram[lowerBound];
					}
					else{
						if(higherBound<actualNbBins){
							interpolatedHistogram[i+translationBin] = (1-Math.abs((float)i-(float)lowerBound*(float)nbBins/(float)actualNbBins)/((float)nbBins/(float)actualNbBins))*(histogram[lowerBound])
									+ (1-Math.abs((float)i-(float)higherBound*(float)nbBins/(float)actualNbBins) / ((float)nbBins/(float)actualNbBins)) *(histogram[higherBound]);
						}
						else{
							interpolatedHistogram[i+translationBin] = histogram[actualNbBins-1];
						}
					}
				}
				else{
					interpolatedHistogram[i] = histogram[actualNbBins-1];
				}
			}
		}
		else{
			outputHistogram = histogram;
			interpolatedHistogram = histogram;
		}

		return outputHistogram;
	}

	// compute Gaussian kernel for density estimation
	float[] computeGaussianKernel(int nbBins,float bandwidth,float maxValue){

		float[] GaussianKernel = new float[nbBins];
		for(int i=0;i<nbBins;i++){
			float correspondingValue=(((float)(i)-(float)(nbBins)/2)*maxValue/(float)(nbBins));
			GaussianKernel[i] = (float)(Math.exp(-Math.pow(correspondingValue/bandwidth,2.)/2.)/(bandwidth*Math.sqrt(2*Math.PI)));
		}

		return GaussianKernel;
	}


	// compute density
	float[] computeDensity(int nbBins,float[] histogram,float[] GaussianKernel){
		float[] density = new float[nbBins];
		for(int i=0;i<nbBins;i++){
			float GaussianNormalization=0;
			for(int j=0;j<nbBins;j++){
				int currentIndex = i+j-nbBins/2;
				if((currentIndex>=0)&&(currentIndex<nbBins)){
					density[i] += histogram[currentIndex]*GaussianKernel[j];
					GaussianNormalization += GaussianKernel[j];
				}
			}
			density[i] /= GaussianNormalization;
		}

		return density;
	}

	// compute histogram for radius
	float[] computeSymmetricalHistogram(float[] input,float[] weight,float[] length,int nbBins,float minValue,float maxValue,float[] interpolatedHistogram){

		int actualNbBins;
		if(nbBins>(2*maxValue+1.)){actualNbBins = (int)Math.ceil(2*maxValue+1);}
		else{actualNbBins = nbBins;}
		
		float[] histogram = new float[actualNbBins],
				outputHistogram = new float[nbBins];
		int[] nbEltsPerBin = new int[nbBins];
		for(int i=0;i<input.length;i++){
			if(!Float.isNaN(input[i])){
				int currentBin = Math.round((float)(actualNbBins-1)*(input[i]-minValue)/(2*maxValue));
				if(currentBin==actualNbBins){currentBin=0;}
				if(!Float.isNaN(weight[i])&&!Float.isNaN(length[i])&&(weight[i]<1000000)){
					histogram[currentBin] += weight[i];
					nbEltsPerBin[i]++;
				}
			}
		}
		
		double totalSum=0.;
		for(int i=0;i<actualNbBins;i++){
			if(nbEltsPerBin[i]>0){
				histogram[i] /= (float)nbEltsPerBin[i];
			}
		}

		if(actualNbBins!=nbBins){
			int translationBin = (int)((float)nbBins/(2*(float)actualNbBins));
			for(int i=0;i<translationBin;i++){
				interpolatedHistogram[i] = histogram[0];
			}
			for(int i=0;i<nbBins;i++){
				outputHistogram[i] = histogram[(int)((float)i*(float)actualNbBins/(float)nbBins)];
				int lowerBound=(int)Math.floor((float)i*(float)actualNbBins/(float)nbBins),
						higherBound=(int)Math.ceil((float)i*(float)actualNbBins/(float)nbBins);
				if(i<(nbBins-translationBin)){
					if(lowerBound==higherBound){
						interpolatedHistogram[i+translationBin] = histogram[lowerBound];
					}
					else{
						if(higherBound<actualNbBins){
							interpolatedHistogram[i+translationBin] = (1-Math.abs((float)i-(float)lowerBound*(float)nbBins/(float)actualNbBins)/((float)nbBins/(float)actualNbBins))*(histogram[lowerBound])
									+ (1-Math.abs((float)i-(float)higherBound*(float)nbBins/(float)actualNbBins) / ((float)nbBins/(float)actualNbBins)) *(histogram[higherBound]);
						}
						else{
							interpolatedHistogram[i+translationBin] = histogram[actualNbBins-1];
						}
					}
				}
				else{
					interpolatedHistogram[i] = histogram[actualNbBins-1];
				}
			}
		}
		else{
			outputHistogram = histogram;
			interpolatedHistogram = histogram;
		}

		return outputHistogram;
	}

	// compute Gaussian kernel for density estimation
	float[] computeSymmetricalGaussianKernel(int nbBins,float bandwidth,float maxValue){

		float[] GaussianKernel = new float[nbBins];
		for(int i=0;i<nbBins;i++){
			float correspondingValue=(((float)(i)-(float)(nbBins)/2)*(2*maxValue+1)/(float)(nbBins));
			GaussianKernel[i] = (float)(Math.exp(-Math.pow(correspondingValue/bandwidth,2.)/2.)/(bandwidth*Math.sqrt(2*Math.PI)));
		}

		return GaussianKernel;
	}

		
	
	// Modified Bessel function of order 0
	double bessel_i0(double X) {
		double
		P1=1.0, P2=3.5156229, P3=3.0899424, P4=1.2067492,
		P5=0.2659732, P6=0.360768e-1, P7=0.45813e-2,
		Q1=0.39894228, Q2=0.1328592e-1, Q3=0.225319e-2,
		Q4=-0.157565e-2, Q5=0.916281e-2, Q6=-0.2057706e-1,
		Q7=0.2635537e-1, Q8=-0.1647633e-1, Q9=0.392377e-2,
		AX = Math.abs(X);
		if (AX < 3.75) {
			double Y = (X/3.75)*(X/3.75);
			return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
		}
		else {
			double Y = 3.75 / AX, BX = Math.exp(AX)/Math.sqrt(AX);
			return (BX * (Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))));
		}
	}

	// Modified Bessel function of order 1
	double bessel_i1(double X) {
		double
		P1=0.5, P2=0.87890594, P3=0.51498869, P4=0.15084934,
		P5=0.2658733e-1, P6=0.301532e-2, P7=0.32411e-3,
		Q1=0.39894228, Q2=-0.3988024e-1, Q3=-0.362018e-2,
		Q4=0.163801e-2, Q5=-0.1031555e-1, Q6=0.2282967e-1,
		Q7=-0.2895312e-1, Q8=0.1787654e-1, Q9=-0.420059e-2,
		AX = Math.abs(X);
		if (AX < 3.75) {
			double Y = (X / 3.75)*(X / 3.75);
			return (X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
		} else {
			double Y = 3.75 / AX, BX = Math.exp(AX)/Math.sqrt(AX);
			return (BX * (Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))));
		}
	}

	// Modified Bessel function
	double bessel_i(int N, double X) {
		double BIGNO = 1e10, BIGNI = 1e-10;
		int IACC = 40, M = (int) (2*((N+Math.floor(Math.sqrt(IACC*N)))));
		if (N==0)  return (bessel_i0(X));
		if (N==1)  return (bessel_i1(X));
		if (X==0.0) return 0.0;
		double TOX = 2.0/X, BIP = 0.0, BI  = 1.0, BSI = 0.0;
		for (int J = M; J>0; J--) {
			double BIM = BIP+J*TOX*BI;
			BIP = BI;
			BI  = BIM;
			if (Math.abs(BI) > BIGNO) {
				BI  = BI*BIGNI;
				BIP = BIP*BIGNI;
				BSI = BSI*BIGNI;
			}
			if (J==N)  BSI = BIP;
		}
		return (BSI*bessel_i0(X)/BI);
	}

	
	//! Rule of Thumb estimation of the concentration for von Mises distribution kernel
	/**
	   \return the concentration for the van mises distribution
	   \param S a vector of angles
	 **/
	float compute_rule_of_thumb_concentration(float[] input) {

		float mu;

		// Compute R
		double kappaRef=0.;
		for(int p=1;p<=3;p++){
			double sc = 0, ss = 0, n = (double)input.length, R=0.;
			for(int i=0;i<input.length;i++){
				double theta = input[i];
				sc += Math.cos(p*theta);
				ss += Math.sin(p*theta);
			}
			double m = Math.atan2(ss,sc);
			for(int i=0;i<input.length;i++){
				double theta = input[i];
				R += Math.cos(p*theta-m);
			}
			R /= n;
			// Estimate the concentration kappa of the data using a newton method
			// http://en.wikipedia.org/wiki/Von_Mises%E2%80%93Fisher_distribution
			double kappa = R * (2.0 - R*R) / (1.0 - R*R), dkappa = 1.0;
			for (int i = 0; i < 1000 && Math.abs(dkappa) > 1e-12; i++) {
				double A = bessel_i(p,kappa) / bessel_i0(kappa);
				dkappa = (A-R) / (1.0 - A * A - A / kappa);
				kappa -= dkappa;
			}
			if(kappa>kappaRef){kappaRef = kappa;}
		}
		// from (Taylor, 2008) deduce the plugin concentration for smoothing
		mu = (float)Math.pow((3.0*(double)(input.length)*kappaRef*kappaRef*bessel_i(2,2.0*kappaRef)) / (4.0*Math.sqrt(Math.PI) * bessel_i0(kappaRef)*bessel_i0(kappaRef)), 2.0/5.0);
		if(Float.isNaN(mu)){mu=0;}

		return mu;

	}


	// compute histogram for theta distribution over several experiments
	float[] computeCircularHistogram(float[] angle,float[] weight,float[] length,int nbBins){

		float[] angleHistogram = new float[nbBins];
		int[] nbEltsPerBin = new int[nbBins];
		for(int i=0;i<angle.length;i++){
			if(!Float.isNaN(angle[i])){
	            int currentBin = (int)Math.round((float)nbBins*angle[i]/(2*Math.PI));
	            if(currentBin==nbBins){currentBin=0;}
	            if(!Float.isNaN(weight[i])&&!Float.isNaN(length[i])){
	            	angleHistogram[currentBin] += weight[i];
	            	nbEltsPerBin[currentBin]++;
	            }
			}
		}

		for(int i=0;i<nbBins;i++){
			if(nbEltsPerBin[i]>0){
				angleHistogram[i] /= (float)nbEltsPerBin[i];
			}
		}
		return angleHistogram;
	}

	// compute von Mises kernel for circular density estimation
	float[] computeVonMisesKernel(int nbBins,float concentration){

		float[] vonMisesKernel = new float[nbBins];
		if(concentration>709.){concentration = 709;}
		for(int i=0;i<nbBins;i++){
			float correspondingValue=(((float)i-(float)nbBins/2)*2*(float)Math.PI/(float)nbBins);
			vonMisesKernel[i] = (float)(Math.exp(concentration*Math.cos(correspondingValue))/(2*Math.PI*bessel_i0(concentration)));
		}

		return vonMisesKernel;
	}

	// compute circular density
	float[] computeCircularDensity(int nbBins,float[] histogram,float[] vonMisesKernel){

		float[] density = new float[nbBins];
		float totalDensitySum=0;
		for(int i=0;i<nbBins;i++){
			for(int j=0;j<nbBins;j++){
				int currentIndex = i+j-nbBins/2;
				if(currentIndex<0){currentIndex += nbBins;}
				if(currentIndex>=nbBins){currentIndex -= nbBins;}
				density[i] += histogram[currentIndex]*vonMisesKernel[j];
			}
			totalDensitySum += density[i];
		}
		for(int i=0;i<nbBins;i++){
			density[i] = density[i]/totalDensitySum;
		}

		return density;
	}

	// export histograms and densities as xls file
	void exportCylindricalFiles(float[] histogram1,float[] histogram2,float[] histogram3,float[] density1,float[] density2,float[] density3,File f,Sequence seq){
		if(f.exists()){
			WritableWorkbook wb = null;
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet depthHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			WritableSheet depthDensitySheet = null;
			int radiusRows=0,angleRows=0,depthRows=0,currentColumn=0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite(f);
				if(wb.getNumberOfSheets()!=6){MessageDialog.showDialog("The excel output file should contain 6 worksheets: Radius histogram, Polar angle histogram, Depth histogram, Radius density, Polar angle density and Depth density.");
				try
				{
					XLSUtil.saveAndClose(wb);
				}
				catch (Exception e)
				{
					throw new IcyHandledException(e.getMessage());
				}
				return;}
				radiusHistogramSheet = wb.getSheet(0);
				angleHistogramSheet = wb.getSheet(1);
				depthHistogramSheet = wb.getSheet(2);
				radiusDensitySheet = wb.getSheet(3);
				angleDensitySheet = wb.getSheet(4);
				depthDensitySheet = wb.getSheet(5);
				currentColumn = radiusHistogramSheet.getColumns();
				if((angleHistogramSheet.getColumns()!=currentColumn)||(depthHistogramSheet.getColumns()!=currentColumn)||(radiusDensitySheet.getColumns()!=currentColumn)||(angleDensitySheet.getColumns()!=currentColumn)||(depthDensitySheet.getColumns()!=currentColumn)){
					MessageDialog.showDialog("The excel output file should have the same number of columns in each worksheet.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				radiusRows = radiusHistogramSheet.getRows();
				if(radiusDensitySheet.getRows()!=radiusRows){
					MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				angleRows = angleHistogramSheet.getRows();
				if(angleDensitySheet.getRows()!=angleRows){
					MessageDialog.showDialog("In the excel file, the worksheets Polar angle histogram and Polar angle density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				depthRows = depthHistogramSheet.getRows();
				if(depthDensitySheet.getRows()!=depthRows){
					MessageDialog.showDialog("In the excel file, the worksheets Depth histogram and Depth density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}

			if(histogram1.length!=(radiusRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(radiusHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(radiusDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, currentColumn, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, currentColumn, bin+1, density1[bin]);
			}
			if(histogram2.length!=(angleRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Polar angle histogram and Polar angle density do not have the same number of rows than the number of bins for polar angle.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(angleHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(angleDensitySheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, 1, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, currentColumn, 1, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, bin+1, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, currentColumn, bin+1, density2[histogram2.length-bin]);
			}
			if(histogram3.length!=(depthRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Depth histogram and Depth density do not have the same number of rows than the number of bins for depth.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(depthHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(depthDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(depthHistogramSheet, currentColumn, bin+1, histogram3[bin]);
				XLSUtil.setCellNumber(depthDensitySheet, currentColumn, bin+1, density3[bin]);
			}
			
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
		else{
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet depthHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			WritableSheet depthDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook(f);
				radiusHistogramSheet = XLSUtil.createNewPage(wb, "Radius histogram");
				angleHistogramSheet = XLSUtil.createNewPage(wb, "Polar angle histogram");
				depthHistogramSheet = XLSUtil.createNewPage(wb, "Depth histogram");
				radiusDensitySheet = XLSUtil.createNewPage(wb, "Radius density");
				angleDensitySheet = XLSUtil.createNewPage(wb, "Polar angle density");
				depthDensitySheet = XLSUtil.createNewPage(wb, "Depth density");
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			XLSUtil.setCellString(radiusHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(radiusDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, 0, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, 0, bin+1, density1[bin]);
			}
			XLSUtil.setCellString(angleHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(angleDensitySheet, 0, 0, seq.getName());
			XLSUtil.setCellNumber(angleHistogramSheet, 0, 1, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, 0, 1, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, 0, bin+1, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, 0, bin+1, density2[histogram2.length-bin]);
			}
			XLSUtil.setCellString(depthHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(depthDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(depthHistogramSheet, 0, bin+1, histogram3[bin]);
				XLSUtil.setCellNumber(depthDensitySheet, 0, bin+1, density3[bin]);
			}
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
	}

	// export histogram or density as xls file
	void exportCylindricalFiles(float[] histogram1,float[] histogram2,float[] density1,float[] density2,File f,Sequence seq){
		if(f.exists()){
			WritableWorkbook wb = null;
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			int radiusRows = 0,angleRows=0,currentColumn=0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite(f);
				if(wb.getNumberOfSheets()!=4){MessageDialog.showDialog("The excel output file should contain 4 worksheets: Radius histogram, Polar angle histogram, Radius density and Polar angle density.");
				try
				{
					XLSUtil.saveAndClose(wb);
				}
				catch (Exception e)
				{
					throw new IcyHandledException(e.getMessage());
				}
				return;}
				radiusHistogramSheet = wb.getSheet(0);
				angleHistogramSheet = wb.getSheet(1);
				radiusDensitySheet = wb.getSheet(2);
				angleDensitySheet = wb.getSheet(3);
				currentColumn = radiusHistogramSheet.getColumns();
				if((angleHistogramSheet.getColumns()!=currentColumn)||(radiusDensitySheet.getColumns()!=currentColumn)||(angleDensitySheet.getColumns()!=currentColumn)){
					MessageDialog.showDialog("The excel output file should have the same number of columns in each worksheet.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				radiusRows = radiusHistogramSheet.getRows();
				if(radiusDensitySheet.getRows()!=radiusRows){
					MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				angleRows = angleHistogramSheet.getRows();
				if(angleDensitySheet.getRows()!=angleRows){
					MessageDialog.showDialog("In the excel file, the worksheets Polar angle histogram and Polar angle density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}

			if(histogram1.length!=(radiusRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(radiusHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(radiusDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, currentColumn, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, currentColumn, bin+1, density1[bin]);
			}
			if(histogram2.length!=(angleRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Polar angle histogram and Polar angle density do not have the same number of rows than the number of bins for polar angle.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(angleHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(angleDensitySheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, 1, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, currentColumn, 1, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, bin+1, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, currentColumn, bin+1, density2[histogram2.length-bin]);
			}

			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
		else{
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook(f);
				radiusHistogramSheet = XLSUtil.createNewPage(wb, "Radius histogram");
				angleHistogramSheet = XLSUtil.createNewPage(wb, "Polar angle histogram");
				radiusDensitySheet = XLSUtil.createNewPage(wb, "Radius density");
				angleDensitySheet = XLSUtil.createNewPage(wb, "Polar angle density");
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			XLSUtil.setCellString(radiusHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(radiusDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, 0, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, 0, bin+1, density1[bin]);
			}
			XLSUtil.setCellString(angleHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(angleDensitySheet, 0, 0, seq.getName());
			XLSUtil.setCellNumber(angleHistogramSheet, 0, 1, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, 0, 1, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, 0, bin+1, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, 0, bin+1, density2[histogram2.length-bin]);
			}

			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
	}

	// export histograms and densities as xls file
	void exportCartesianFiles(float[] histogram1,float[] histogram2,float[] histogram3,float[] density1,float[] density2,float[] density3,File f,Sequence seq){
		if(f.exists()){
			WritableWorkbook wb = null;
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet zHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			WritableSheet zDensitySheet = null;
			int xRows=0,yRows=0,zRows=0,currentColumn=0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite(f);
				if(wb.getNumberOfSheets()!=6){MessageDialog.showDialog("The excel output file should contain 6 worksheets: X histogram, Y histogram, Z histogram, X density, Y density and Z density.");
				try
				{
					XLSUtil.saveAndClose(wb);
				}
				catch (Exception e)
				{
					throw new IcyHandledException(e.getMessage());
				}
				return;}
				xHistogramSheet = wb.getSheet(0);
				yHistogramSheet = wb.getSheet(1);
				zHistogramSheet = wb.getSheet(2);
				xDensitySheet = wb.getSheet(3);
				yDensitySheet = wb.getSheet(4);
				zDensitySheet = wb.getSheet(5);
				currentColumn = xHistogramSheet.getColumns();
				if((yHistogramSheet.getColumns()!=currentColumn)||(zHistogramSheet.getColumns()!=currentColumn)||(xDensitySheet.getColumns()!=currentColumn)||(yDensitySheet.getColumns()!=currentColumn)||(zDensitySheet.getColumns()!=currentColumn)){
					MessageDialog.showDialog("The excel output file should have the same number of columns in each worksheet.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				xRows = xHistogramSheet.getRows();
				if(xDensitySheet.getRows()!=xRows){
					MessageDialog.showDialog("In the excel file, the worksheets X histogram and X density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				yRows = yHistogramSheet.getRows();
				if(yDensitySheet.getRows()!=yRows){
					MessageDialog.showDialog("In the excel file, the worksheets Y histogram and Y density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				zRows = zHistogramSheet.getRows();
				if(zDensitySheet.getRows()!=zRows){
					MessageDialog.showDialog("In the excel file, the worksheets Z histogram and Z density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}

			if(histogram1.length!=(xRows-1)){MessageDialog.showDialog("In the excel file, the worksheets X histogram and X density do not have the same number of rows than the number of bins for abscissa.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(xHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(xDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, currentColumn, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, currentColumn, bin+1, density1[bin]);
			}
			if(histogram2.length!=(yRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Y histogram and Y density do not have the same number of rows than the number of bins for ordinate.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(yHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(yDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, currentColumn, bin+1, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, currentColumn, bin+1, density2[bin]);
			}
			if(histogram3.length!=(zRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Z histogram and Z density do not have the same number of rows than the number of bins for height.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(zHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(zDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(zHistogramSheet, currentColumn, bin+1, histogram3[bin]);
				XLSUtil.setCellNumber(zDensitySheet, currentColumn, bin+1, density3[bin]);
			}
			
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
		else{
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet zHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			WritableSheet zDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook(f);
				xHistogramSheet = XLSUtil.createNewPage(wb, "X histogram");
				yHistogramSheet = XLSUtil.createNewPage(wb, "Y histogram");
				zHistogramSheet = XLSUtil.createNewPage(wb, "Z histogram");
				xDensitySheet = XLSUtil.createNewPage(wb, "X density");
				yDensitySheet = XLSUtil.createNewPage(wb, "Y density");
				zDensitySheet = XLSUtil.createNewPage(wb, "Z density");
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			XLSUtil.setCellString(xHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(xDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, 0, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, 0, bin+1, density1[bin]);
			}
			XLSUtil.setCellString(yHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(yDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, 0, bin+1, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, 0, bin+1, density2[bin]);
			}
			XLSUtil.setCellString(zHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(zDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(zHistogramSheet, 0, bin+1, histogram3[bin]);
				XLSUtil.setCellNumber(zDensitySheet, 0, bin+1, density3[bin]);
			}
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
	}

	// export histogram or density as xls file
	void exportCartesianFiles(float[] histogram1,float[] histogram2,float[] density1,float[] density2,File f,Sequence seq){
		if(f.exists()){
			WritableWorkbook wb = null;
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			int xRows=0,yRows=0,currentColumn=0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite(f);
				if(wb.getNumberOfSheets()!=4){MessageDialog.showDialog("The excel output file should contain 4 worksheets: X histogram, Y histogram, X density and Y density.");
				try
				{
					XLSUtil.saveAndClose(wb);
				}
				catch (Exception e)
				{
					throw new IcyHandledException(e.getMessage());
				}
				return;}
				xHistogramSheet = wb.getSheet(0);
				yHistogramSheet = wb.getSheet(1);
				xDensitySheet = wb.getSheet(2);
				yDensitySheet = wb.getSheet(3);
				currentColumn = xHistogramSheet.getColumns();
				if((yHistogramSheet.getColumns()!=currentColumn)||(xDensitySheet.getColumns()!=currentColumn)||(yDensitySheet.getColumns()!=currentColumn)){
					MessageDialog.showDialog("The excel output file should have the same number of columns in each worksheet.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				xRows = xHistogramSheet.getRows();
				if(xDensitySheet.getRows()!=xRows){
					MessageDialog.showDialog("In the excel file, the worksheets X histogram and X density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				yRows = yHistogramSheet.getRows();
				if(yDensitySheet.getRows()!=yRows){
					MessageDialog.showDialog("In the excel file, the worksheets Y histogram and Y density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}

			if(histogram1.length!=(xRows-1)){MessageDialog.showDialog("In the excel file, the worksheets X histogram and X density do not have the same number of rows than the number of bins for abscissa.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(xHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(xDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, currentColumn, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, currentColumn, bin+1, density1[bin]);
			}
			if(histogram2.length!=(yRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Y histogram and Y density do not have the same number of rows than the number of bins for ordinate.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(yHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(yDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, currentColumn, bin+1, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, currentColumn, bin+1, density2[bin]);
			}
			
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
		else{
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook(f);
				xHistogramSheet = XLSUtil.createNewPage(wb, "X histogram");
				yHistogramSheet = XLSUtil.createNewPage(wb, "Y histogram");
				xDensitySheet = XLSUtil.createNewPage(wb, "X density");
				yDensitySheet = XLSUtil.createNewPage(wb, "Y density");
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			XLSUtil.setCellString(xHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(xDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, 0, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, 0, bin+1, density1[bin]);
			}
			XLSUtil.setCellString(yHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(yDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, 0, bin+1, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, 0, bin+1, density2[bin]);
			}
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
	}

	// export histograms and densities as xls file
	void exportSphericalFiles(float[] histogram1,float[] histogram2,float[] histogram3,float[] density1,float[] density2,float[] density3,File f,Sequence seq){
		if(f.exists()){
			WritableWorkbook wb = null;
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angle1HistogramSheet = null;
			WritableSheet angle2HistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angle1DensitySheet = null;
			WritableSheet angle2DensitySheet = null;
			int radiusRows=0,angle1Rows=0,angle2Rows=0,currentColumn=0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite(f);
				if(wb.getNumberOfSheets()!=6){MessageDialog.showDialog("The excel output file should contain 6 worksheets: Radius histogram, Colatitude histogram, Azimuth angle histogram, Radius density, Colatitude density and Azimuth angle density.");
				try
				{
					XLSUtil.saveAndClose(wb);
				}
				catch (Exception e)
				{
					throw new IcyHandledException(e.getMessage());
				}
				return;}
				radiusHistogramSheet = wb.getSheet(0);
				angle1HistogramSheet = wb.getSheet(1);
				angle2HistogramSheet = wb.getSheet(2);
				radiusDensitySheet = wb.getSheet(3);
				angle1DensitySheet = wb.getSheet(4);
				angle2DensitySheet = wb.getSheet(5);
				currentColumn = radiusHistogramSheet.getColumns();
				if((angle1HistogramSheet.getColumns()!=currentColumn)||(angle2HistogramSheet.getColumns()!=currentColumn)||(radiusDensitySheet.getColumns()!=currentColumn)||(angle2DensitySheet.getColumns()!=currentColumn)||(angle2DensitySheet.getColumns()!=currentColumn)){
					MessageDialog.showDialog("The excel output file should have the same number of columns in each worksheet.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				radiusRows = radiusHistogramSheet.getRows();
				if(radiusDensitySheet.getRows()!=radiusRows){
					MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				angle1Rows = angle1HistogramSheet.getRows();
				if(angle1DensitySheet.getRows()!=angle1Rows){
					MessageDialog.showDialog("In the excel file, the worksheets Colatitude histogram and Colatitude density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				angle2Rows = angle2HistogramSheet.getRows();
				if(angle2DensitySheet.getRows()!=angle2Rows){
					MessageDialog.showDialog("In the excel file, the worksheets Azimuth angle histogram and Azimuth angle density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}

			if(histogram1.length!=(radiusRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(radiusHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(radiusDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, currentColumn, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, currentColumn, bin+1, density1[bin]);
			}
			if(histogram2.length!=(angle1Rows-1)){MessageDialog.showDialog("In the excel file, the worksheets Colatitude histogram and Colatitude density do not have the same number of rows than the number of bins for colatitude.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(angle1HistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(angle1DensitySheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellNumber(angle1HistogramSheet, currentColumn, 1, histogram2[0]);
			XLSUtil.setCellNumber(angle1DensitySheet, currentColumn, 1, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angle1HistogramSheet, currentColumn, bin+1, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angle1DensitySheet, currentColumn, bin+1, density2[histogram2.length-bin]);
			}
			if(histogram3.length!=(angle2Rows-1)){MessageDialog.showDialog("In the excel file, the worksheets Azimuth angle histogram and Azimuth angle density do not have the same number of rows than the number of bins for azimuth angle.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(angle2HistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(angle2DensitySheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellNumber(angle2HistogramSheet, currentColumn, 1, histogram3[0]);
			XLSUtil.setCellNumber(angle2DensitySheet, currentColumn, 1, density3[0]);
			for(int bin=1;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(angle2HistogramSheet, currentColumn, bin+1, histogram2[histogram3.length-bin]);
				XLSUtil.setCellNumber(angle2DensitySheet, currentColumn, bin+1, density2[histogram3.length-bin]);
			}

			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
		else{
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angle1HistogramSheet = null;
			WritableSheet angle2HistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angle1DensitySheet = null;
			WritableSheet angle2DensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook(f);
				radiusHistogramSheet = XLSUtil.createNewPage(wb, "Radius histogram");
				angle1HistogramSheet = XLSUtil.createNewPage(wb, "Azimuth angle histogram");
				angle2HistogramSheet = XLSUtil.createNewPage(wb, "Colatitude histogram");
				radiusDensitySheet = XLSUtil.createNewPage(wb, "Radius density");
				angle1DensitySheet = XLSUtil.createNewPage(wb, "Azimuth angle density");
				angle2DensitySheet = XLSUtil.createNewPage(wb, "Colatitude density");
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			XLSUtil.setCellString(radiusHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(radiusDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, 0, bin+1, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, 0, bin+1, density1[bin]);
			}
			XLSUtil.setCellString(angle1HistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(angle1DensitySheet, 0, 0, seq.getName());
			XLSUtil.setCellNumber(angle1HistogramSheet, 0, 1, histogram2[0]);
			XLSUtil.setCellNumber(angle1DensitySheet, 0, 1, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angle1HistogramSheet, 0, bin+1, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angle1DensitySheet, 0, bin+1, density2[histogram2.length-bin]);
			}
			XLSUtil.setCellString(angle2HistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(angle2DensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(angle2HistogramSheet, 0, bin+1, histogram3[bin]);
				XLSUtil.setCellNumber(angle2DensitySheet, 0, bin+1, density3[bin]);
			}
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
	}

	// export histogram or density as xls file
	void exportDistanceFiles(float[] histogram,float[] density,File f,Sequence seq){
		if(f.exists()){
			WritableWorkbook wb = null;
			WritableSheet distanceHistogramSheet = null;
			WritableSheet distanceDensitySheet = null;
			int distanceRows = 0,currentColumn=0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite(f);
				if(wb.getNumberOfSheets()!=2){MessageDialog.showDialog("The excel output file should contain 2 worksheets: Distance to cell border histogram and Distance to cell border density.");
				try
				{
					XLSUtil.saveAndClose(wb);
				}
				catch (Exception e)
				{
					throw new IcyHandledException(e.getMessage());
				}
				return;}
				distanceHistogramSheet = wb.getSheet(0);
				distanceDensitySheet = wb.getSheet(1);
				currentColumn = distanceHistogramSheet.getColumns();
				if(distanceDensitySheet.getColumns()!=currentColumn){
					MessageDialog.showDialog("The excel output file should have the same number of columns in each worksheet.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
				distanceRows = distanceHistogramSheet.getRows();
				if(distanceDensitySheet.getRows()!=distanceRows){
					MessageDialog.showDialog("In the excel file, the worksheets Distance to cell border histogram and Distance to cell border density should have the same number of rows.");
					try
					{
						XLSUtil.saveAndClose(wb);
					}
					catch (Exception e)
					{
						throw new IcyHandledException(e.getMessage());
					}
					return;}
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}

			if(histogram.length!=(distanceRows-1)){MessageDialog.showDialog("In the excel file, the worksheets Distance to cell border histogram and Distance to cell border density do not have the same number of rows than the number of bins for distance to cell border.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellString(distanceHistogramSheet, currentColumn, 0, seq.getName());
			XLSUtil.setCellString(distanceDensitySheet, currentColumn, 0, seq.getName());
			for(int bin=0;bin<histogram.length;bin++){
				XLSUtil.setCellNumber(distanceHistogramSheet, currentColumn, bin+1, histogram[bin]);
				XLSUtil.setCellNumber(distanceDensitySheet, currentColumn, bin+1, density[bin]);
			}

			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
		else{
			WritableSheet distanceHistogramSheet = null;
			WritableSheet distanceDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook(f);
				distanceHistogramSheet = XLSUtil.createNewPage(wb, "Distance histogram");
				distanceDensitySheet = XLSUtil.createNewPage(wb, "Distance density");
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			XLSUtil.setCellString(distanceHistogramSheet, 0, 0, seq.getName());
			XLSUtil.setCellString(distanceDensitySheet, 0, 0, seq.getName());
			for(int bin=0;bin<histogram.length;bin++){
				XLSUtil.setCellNumber(distanceHistogramSheet, 0, bin+1, histogram[bin]);
				XLSUtil.setCellNumber(distanceDensitySheet, 0, bin+1, density[bin]);
			}

			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
		}
	}

	// compute CEMD
	float[][] CEMD(float[][] density,int nbExperiments,int nbBins){
		float[][] output = new float[nbExperiments][nbExperiments];
		float[][][] intermediateDensities = new float[nbExperiments][nbBins][nbBins];
		for(int xp=0;xp<nbExperiments;xp++){
			for(int k=0;k<nbBins;k++){
				intermediateDensities[xp][k][0] = density[xp][k];
					for(int bin=1;bin<nbBins;bin++){
						int currentBin = bin+k;
						if(currentBin>=nbBins){currentBin -= nbBins;}
						intermediateDensities[xp][k][bin] = intermediateDensities[xp][k][bin-1] + density[xp][currentBin];
				}
			}
		}

		for(int i=0;i<nbExperiments;i++){
			for(int u=i+1;u<nbExperiments;u++){
				float minOutput=0;
				for(int bin=0;bin<nbBins;bin++){
					minOutput += Math.abs(intermediateDensities[i][0][bin]-intermediateDensities[u][0][bin]);
				}
				for(int k=1;k<nbBins;k++){
					float currentOutput=0;
					for(int bin=0;bin<nbBins;bin++){
						currentOutput += Math.abs(intermediateDensities[i][k][bin]-intermediateDensities[u][k][bin]);
					}
					if(currentOutput<minOutput){minOutput = currentOutput;}
				}
				output[i][u] = minOutput;
				output[u][i] = output[i][u];
			}
		}
		return output;
	}


	// compute EMD
	float[][] EMD(float[][] density,int nbExperiments,int nbBins){
		float[][] output = new float[nbExperiments][nbExperiments],cumulatedDensities = density;
		for(int bin=1;bin<nbBins;bin++){
			for(int xp=0;xp<nbExperiments;xp++){
				cumulatedDensities[xp][bin] += cumulatedDensities[xp][bin-1];
			}
		}
		for(int i=0;i<nbExperiments;i++){
			for(int u=i+1;u<nbExperiments;u++){
				for(int bin=0;bin<nbBins;bin++){
					output[i][u] += Math.abs(cumulatedDensities[i][bin]-cumulatedDensities[u][bin]);
				}
				output[u][i] = output[i][u];
			}
		}
		return output;
	}

	
	@Override
	protected void initialize()
	{
		// number of bins
		// Cartesian coordinate system
		int nbBinsForX = 100;
		int nbBinsForY = 100;
		int nbBinsForCartesianDepth = 100;
		// cylindrical coordinate system
		int nbBinsForCylindricalRadius = 100;
		int nbBinsForAngle = 180;
		int nbBinsForCylindricalDepth= 100;
		// spherical coordinate system
		int nbBinsForSphericalRadius = 100;
		int nbBinsForFirstAngle = 180;
		int nbBinsForSecondAngle = 180;
		// distance distribution (distance to cell border for example)
		int nbBinsForDistance = 100;
		
		// coordinate system center
		// Cartesian coordinate system
		int ref1Xcartesian = -1;
		int ref1Ycartesian = -1;
		// cylindrical coordinate system
		int ref1Xcylindrical = -1;
		int ref1Ycylindrical = -1;
		// spherical coordinate system
		int ref1Xspherical = -1;
		int ref1Yspherical= -1;
		int ref1Zspherical = -1;
		
		// second point to define coordinate system reference direction
		// Cartesian coordinate system
		int ref2Xcartesian = -1;
		int ref2Ycartesian = -1;
		// cylindrical coordinate system
		int ref2Xcylindrical = -1;
		int ref2Ycylindrical = -1;
		// spherical coordinate system
		int ref2Xspherical = -1;
		int ref2Yspherical = -1;

		// coordinate system
		String[] coordinateSystemPossibilities = {"Cylindrical","Spherical","Cartesian","Distance to cell border"};
		String coordinateSystem = coordinateSystemPossibilities[0];
		
		// input data
		input;
		
		// input cell mask for cell normalization (and forbidden region if defined)
		image cellMaskDensities = new image;
		image forbiddenRegionDensities = new image;

		// output xls files
		File exportCylindricalExcelFile = new File("Output cylindrical xls file", "");
		File exportCartesianExcelFile = new File("Output Cartesian xls file", "");
		File exportSphericalExcelFile = new File("Output spherical xls file", "");
		File exportDistanceExcelFile = new File("Output distance to cell border xls file", "");
		
	}
	
	@Override
	public void execute() {

		final Sequence cellMaskImage = cellMaskDensities.getValue();
		final Sequence forbiddenRegionImage = forbiddenRegionDensities.getValue();

		// test on output xls files
		if(coordinateSystem.getValue()=="Cylindrical"){
			if(exportCylindricalExcelFile.getValue(false)==null){
				MessageDialog.showDialog("You need to specify the output xls file");
				return;
			}
		}
		if(coordinateSystem.getValue()=="Cartesian"){
			if(exportCartesianExcelFile.getValue(false)==null){
				MessageDialog.showDialog("You need to specify the output xls file");
				return;
			}
		}
		if(coordinateSystem.getValue()=="Spherical"){
			if(exportSphericalExcelFile.getValue(false)==null){
				MessageDialog.showDialog("You need to specify the output xls file");
				return;
			}
		}
		if(coordinateSystem.getValue()=="Distance to cell border"){
			if(exportDistanceExcelFile.getValue(false)==null){
				MessageDialog.showDialog("You need to specify the output xls file");
				return;
			}
		}

		// interesting variables
		int arraySize=inputImage.getSizeX()*inputImage.getSizeY(),
				width=inputImage.getSizeX(),
				height=inputImage.getSizeY(),
				depth=inputImage.getSizeZ(),
				nbFrames=inputImage.getSizeT(),
				cellMaskDepth=1,
				forbiddenRegionDepth=1;
		
		// variable initialization
		float[] weight = new float[0],
				component1 = new float[0],
				component2 = new float[0],
				component3 = new float[0],
				distance1 = new float[0],
				distance2 = new float[0];
		int cellCenter1=0,cellCenter2=0,cellCenter3=0,
				referenceCenter1,referenceCenter2,referenceCenter3;
		float referenceDirection1=(float)Math.PI/2;

		// cell mask
		// cell center is computed
		// in case no reference center is given by the user, the cell center is defined as the reference center
		int[][] cellMaskArrayInit = new int[cellMaskDepth][arraySize];
		if(cellMaskImage==null){
			if((coordinateSystem.getValue()=="Cylindrical")||(coordinateSystem.getValue()=="Cartesian")){
				for(int i=0;i<arraySize;i++){
					cellMaskArrayInit[0][i] = 1;
				}
				cellCenter1 = width/2;
				cellCenter2 = height/2;
				cellCenter3 = 0;
			}
			else{
				for(int z=0;z<cellMaskDepth;z++){
					for(int i=0;i<arraySize;i++){
						cellMaskArrayInit[z][i] = 1;
					}
				}
				cellCenter1 = width/2;
				cellCenter2 = height/2;
				cellCenter3 = depth/2;
			}
		}
		else{
			cellMaskArrayInit = Array2DUtil.arrayToIntArray(cellMaskImage.getDataXYZ(0,0), cellMaskImage.isSignedDataType());
			int nbTrajectories=0;
			for(int z=0;z<cellMaskDepth;z++){
				for(int y=0;y<height;y++){
					for(int x=0;x<width;x++){
						if(cellMaskArrayInit[z][y*width+x]>0){
							cellMaskArrayInit[z][y*width+x] = 1;
							cellCenter1 += x;
							cellCenter2 += y;
							cellCenter3 += z;
							nbTrajectories++;
						}
					}
				}
			}
			cellCenter1 = (int)((float)cellCenter1/(float)nbTrajectories);
			cellCenter2 = (int)((float)cellCenter2/(float)nbTrajectories);
			cellCenter3 = (int)((float)cellCenter3/(float)nbTrajectories);
		}
		int[][] cellMaskArray = new int[depth][arraySize];
		if(cellMaskDepth==depth){cellMaskArray = cellMaskArrayInit;}
		else{
			for(int z=0;z<depth;z++){
				for(int i=0;i<arraySize;i++){
					cellMaskArray[z][i] = cellMaskArrayInit[0][i];
				}
			}
		}

		// forbidden region
		int[][] forbiddenRegionArrayInit = new int[forbiddenRegionDepth][arraySize];
		if(forbiddenRegionImage==null){
			if(coordinateSystem.getValue()=="Spherical"){
				for(int z=0;z<forbiddenRegionDepth;z++){
					for(int i=0;i<arraySize;i++){
						forbiddenRegionArrayInit[z][i] = 0;
					}
				}
			}
			else{
				for(int i=0;i<arraySize;i++){
					forbiddenRegionArrayInit[0][i] = 0;
				}
			}
		}
		else{
			forbiddenRegionArrayInit = Array2DUtil.arrayToIntArray(forbiddenRegionImage.getDataXYZ(0,0), forbiddenRegionImage.isSignedDataType());
			for(int z=0;z<forbiddenRegionDepth;z++){
				for(int y=0;y<height;y++){
					for(int x=0;x<width;x++){
						if(forbiddenRegionArrayInit[z][y*width+x]>0){
							forbiddenRegionArrayInit[z][y*width+x] = 1;
						}
					}
				}
			}
		}
		int[][] forbiddenRegionArray = new int[depth][arraySize];
		if(forbiddenRegionDepth==depth){forbiddenRegionArray = forbiddenRegionArrayInit;}
		else{
			for(int z=0;z<depth;z++){
				for(int i=0;i<arraySize;i++){
					forbiddenRegionArray[z][i] = forbiddenRegionArrayInit[0][i];
				}
			}
		}

		// outer circle definition from mask
		int[][] cellBorder = new int[0][0];
		if((coordinateSystem.getValue()=="Cylindrical")||(coordinateSystem.getValue()=="Cartesian")){
			cellBorder = computePseudo3DCellBorder(cellMaskArrayInit,width,height,cellMaskDepth);
		}
		else{
			cellBorder = compute3DCellBorder(cellMaskArrayInit,width,height,cellMaskDepth);
		}

		// inner border if forbidden region and reference center inside forbidden region
		int[][] innerBorder = new int[0][0]; 
		if(forbiddenRegionImage!=null){
			if((coordinateSystem.getValue()=="Cylindrical")||(coordinateSystem.getValue()=="Cartesian")){
				innerBorder = computePseudo3DInnerBorder(forbiddenRegionArrayInit,width,height,forbiddenRegionDepth);
			}
			else{
				innerBorder = compute3DInnerBorder(forbiddenRegionArrayInit,width,height,forbiddenRegionDepth);
			}
		}

		// sum over time and extract intensity
		float[] intensitySumOverTime = new float[arraySize];
		for(int t=0;t<nbFrames;t++){
			float[][] inputArray = Array2DUtil.arrayToFloatArray(inputImage.getDataXYZ(t,channelOfInterestDensities.getValue()), inputImage.isSignedDataType());
			for(int z=0;z<depth;z++){
				for(int y=0;y<height;y++){
					for(int x=0;x<width;x++){
						if((inputArray[z][y*width+x]>0.001)&&(cellMaskArray[z][y*width+x]>0)){
							intensitySumOverTime[y*width+x] += inputArray[z][y*width+x];
						}
					}
				}
			}
		}

		// compute events coordinates
		if(coordinateSystem.getValue()=="Cylindrical"){

			// reference center
			if((ref1Xcylindrical.getValue()>-1)&&(ref1Xcylindrical.getValue()<width)&&(ref1Ycylindrical.getValue()>-1)&&(ref1Ycylindrical.getValue()<height)){
				referenceCenter1 = ref1Xcylindrical.getValue();
				referenceCenter2 = ref1Ycylindrical.getValue();
			}
			else{
				referenceCenter1 = cellCenter1;
				referenceCenter2 = cellCenter2;
			}

			// reference direction
			if((ref2Xcylindrical.getValue()>-1)&&(ref2Xcylindrical.getValue()<width)&&(ref2Ycylindrical.getValue()>-1)&&(ref2Ycylindrical.getValue()<height)){
				referenceDirection1 = (float)(Math.atan2((referenceCenter2-ref2Ycylindrical.getValue()),(ref2Xcylindrical.getValue()-referenceCenter1)) + Math.PI/2.);
			}
			else{
				if((referenceCenter1!=cellCenter1)&&(referenceCenter2!=cellCenter2)){
					referenceDirection1 = (float)(Math.atan2((referenceCenter2-cellCenter2),(cellCenter1-referenceCenter1)) + Math.PI/2.);
				}
			}

			
			float maxRadius=0;
			// number of trajectories 
			int nbTrajectories;
			component1 = new float[nbTrajectories];
			component2 = new float[nbTrajectories];
			component3 = new float[nbTrajectories];
			weight = new float[nbTrajectories];
			distance1 = new float[nbTrajectories];
			distance2 = new float[nbTrajectories];
			distance2 = new float[nbTrajectories];

			// compute normalizing distances for each spatial point so it's done once for all
			float[][] cellSegment = new float[depth][arraySize];
			int[][] insideXCoord = new int[depth][arraySize];
			int[][] insideYCoord = new int[depth][arraySize];
			for(int z=0;z<depth;z++){
				for(int y=0;y<height;y++){
					for(int x=0;x<width;x++){
						if(intensitySumOverTime[y*width+x]>0){
							if((cellMaskImage!=null)||(forbiddenRegionImage!=null)){
								if(forbiddenRegionImage!=null){
									if(forbiddenRegionArray[z][referenceCenter2*width+referenceCenter1]>0){
										cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,innerBorder,width,height,depth,insideXCoord,insideYCoord);
									}
									else{
										cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,width,height,depth);
									}
								}
								else{
									cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,width,height,depth);
								}
							}
							else{
								cellSegment[z][y*width+x] = 1;
							}
						}
					}
				}
			}

			int cpt=0;
			///////////////////////////////////////
			// loop over trajectories
			///////////////////////////////////////
			for(trajectories){
				// extract coordinates x,y,z associated to current trajectory
				// the middle point should be favoured
				if((cellMaskArray[z][y*width+x]>0)&&(forbiddenRegionArray[z][y*width+x]==0)&&(inputArray[z][y*width+x]>0)){
					// compute angle
					float theta=(float)(Math.atan2(referenceCenter2-y,x-referenceCenter1));

					// take into account the orientation reference
					float orientation=(theta-referenceDirection1);
					while(orientation<0){orientation+=(2*Math.PI);}
					while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

					// initialization
					float currentRadius=0;
					if(forbiddenRegionImage!=null){
						if(forbiddenRegionArray[z][referenceCenter2*width+referenceCenter1]>0){
							currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[z][y*width+x]-x,2.)+Math.pow(insideYCoord[z][y*width+x]-y,2.)));
						}
						else{
							currentRadius =	(float)(Math.sqrt(Math.pow(referenceCenter1-x,2.)+Math.pow(referenceCenter2-y,2.)));
						}
					}
					else{
						currentRadius =	(float)(Math.sqrt(Math.pow(referenceCenter1-x,2.)+Math.pow(referenceCenter2-y,2.)));
					}

					// store coordinates 
					if(currentRadius>maxRadius){maxRadius = currentRadius;}
					extractCylindricalCoordinates(inputArray, x, y, z, currentRadius, orientation, cellSegment[z][y*width+x], considerIntensityCylindrical.getValue(), component1, component2, component3, weight, distance1, width, height, depth, cpt);
					cpt++;
			}

			// radius
			// compute density for radius distribution for each experiment
			float bandwidthForGaussianDistributionForRadius = compute_rule_of_thumb_bandwidth(component1);

			// compute radius histograms
			// To modify distanceToZborder -> distanceToPlaneBorder
			float[] radiusHistogram = computeRadiusHistogram(component1,weight,distance1,nbBinsForCylindricalRadius.getValue());

			// compute Gaussian kernel for each experiment
			float[] GaussianKernelForRadius = computeGaussianKernel(nbBinsForCylindricalRadius.getValue(),bandwidthForGaussianDistributionForRadius,maxRadius);

			// compute radius density
			float[] radiusDensity = computeDensity(nbBinsForCylindricalRadius.getValue(),radiusHistogram,GaussianKernelForRadius);

			// theta
			// compute density for theta distribution for each experiment
			float concentrationForVonMisesDistribution = compute_rule_of_thumb_concentration(component2);

			// compute theta histograms
			float[] angleHistogram = computeCircularHistogram(component2,weight,distance1,nbBinsForAngle.getValue());

			// compute von Mises kernel for each experiment
			float[] vonMisesKernel = computeVonMisesKernel(nbBinsForAngle.getValue(),concentrationForVonMisesDistribution);

			// compute circular density
			float[] angleDensity = computeCircularDensity(nbBinsForAngle.getValue(),angleHistogram,vonMisesKernel);

			// depth
			// compute density for depth distribution for each experiment
			if(depth>1){
				// compute density for depth distribution for each experiment
				float bandwidthForGaussianDistributionForDepth = compute_rule_of_thumb_bandwidth(component3);

				// compute depth histograms
				float[] interpolatedDepthHistogram = new float[nbBinsForCylindricalDepth.getValue()],
						depthHistogram = computeDepthHistogram(component3,weight,distance2,nbBinsForCylindricalDepth.getValue(),depth-1,interpolatedDepthHistogram);

				// compute Gaussian kernel for each experiment
				float[] GaussianKernelForDepth = computeGaussianKernel(nbBinsForCylindricalDepth.getValue(),bandwidthForGaussianDistributionForDepth,depth-1);

				// compute depth density
				float[] depthDensity = computeDensity(nbBinsForCylindricalDepth.getValue(),depthHistogram,GaussianKernelForDepth);

				File f = exportCylindricalExcelFile.getValue(true);
				if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
				exportCylindricalFiles(radiusHistogram,angleHistogram,depthHistogram,radiusDensity,angleDensity,depthDensity,f,inputImage);

			}
			else{
				File f = exportCylindricalExcelFile.getValue(true);
				if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
				exportCylindricalFiles(radiusHistogram,angleHistogram,radiusDensity,angleDensity,f,inputImage);
			}
		}

		else{
			// compute events coordinates
			if(coordinateSystem.getValue()=="Cartesian"){

				if((ref1Xcartesian.getValue()>-1)&&(ref1Xcartesian.getValue()<width)&&(ref1Ycartesian.getValue()>-1)&&(ref1Ycartesian.getValue()<height)){
					referenceCenter1 = ref1Xcartesian.getValue();
					referenceCenter2 = ref1Ycartesian.getValue();
				}
				else{
					referenceCenter1 = cellCenter1;
					referenceCenter2 = cellCenter2;
				}

				// reference direction
				if((ref2Xcartesian.getValue()>-1)&&(ref2Xcartesian.getValue()<width)&&(ref2Ycartesian.getValue()>-1)&&(ref2Ycartesian.getValue()<height)){
					referenceDirection1 = (float)(Math.atan2((referenceCenter2-ref2Ycartesian.getValue()),(ref2Xcartesian.getValue()-referenceCenter1)) + Math.PI/2.);
				}
				else{
					if((referenceCenter1!=cellCenter1)&&(referenceCenter2!=cellCenter2)){
						referenceDirection1 = (float)(Math.atan2((referenceCenter2-cellCenter2),(cellCenter1-referenceCenter1)) + Math.PI/2.);
					}
				}

				float minX=1000000,minY=100000,maxX=0,maxY=0;
				int nbTrajectories;
				component1 = new float[nbTrajectories];
				component2 = new float[nbTrajectories];
				component3 = new float[nbTrajectories];
				weight = new float[nbTrajectories];
				distance1 = new float[nbTrajectories];
				distance2 = new float[nbTrajectories];

				// compute normalizing distance
				// compute normalizing distances for each spatial point so it's done once for all
				float[][] cellSegment = new float[depth][arraySize];
				int[][] insideXCoord = new int[depth][arraySize];
				int[][] insideYCoord = new int[depth][arraySize];
				for(int z=0;z<depth;z++){
					for(int y=0;y<height;y++){
						for(int x=0;x<width;x++){
							if(intensitySumOverTime[y*width+x]>0.001){
								if((cellMaskImage!=null)||(forbiddenRegionImage!=null)){
									if(forbiddenRegionImage!=null){
										if(forbiddenRegionArray[z][referenceCenter2*width+referenceCenter1]>0){
											cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,innerBorder,width,height,depth,insideXCoord,insideYCoord);
										}
										else{
											cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,width,height,depth);
										}
									}
									else{
										cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,width,height,depth);
									}
								}
								else{
									cellSegment[z][y*width+x] = 1;
								}
							}
						}
					}
				}

				int cpt=0;
				///////////////////////////////////////
				// loop over trajectories
				///////////////////////////////////////
				for(trajectories){
					// extract coordinates x,y,z associated to current trajectory
					// the middle point should be favoured
					if((cellMaskArray[z][y*width+x]>0)&&(forbiddenRegionArray[z][y*width+x]==0)&&(inputArray[z][y*width+x]>0)){
						// compute angle
						float theta=(float)(Math.atan2(referenceCenter2-y,x-referenceCenter1));
						// take into account the orientation reference
						float orientation=(theta-referenceDirection1);
						while(orientation<0){orientation+=(2*Math.PI);}
						while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

						// initialization
						float currentRadius=0;
						if(forbiddenRegionImage!=null){
							if(forbiddenRegionArray[z][referenceCenter2*width+referenceCenter1]>0){
								currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[z][y*width+x]-x,2.)+Math.pow(insideYCoord[z][y*width+x]-y,2.)));
							}
							else{
								currentRadius =	(float)(Math.sqrt(Math.pow(referenceCenter1-x,2.)+Math.pow(referenceCenter2-y,2.)));
							}
						}
						else{
							currentRadius =	(float)(Math.sqrt(Math.pow(referenceCenter1-x,2.)+Math.pow(referenceCenter2-y,2.)));
						}

						if(Math.abs(currentRadius*(float)Math.cos(orientation))>maxX){maxX = Math.abs(currentRadius*(float)Math.cos(orientation));}
						if(Math.abs(currentRadius*(float)Math.sin(orientation))>maxY){maxY = Math.abs(currentRadius*(float)Math.sin(orientation));}
						if(currentRadius*(float)Math.cos(orientation)<minX){minX = (currentRadius*(float)Math.cos(orientation));}
						if(currentRadius*(float)Math.sin(orientation)<minY){minY = (currentRadius*(float)Math.sin(orientation));}
						extractCartesianCoordinates(inputArray, x, y, z, currentRadius, orientation, cellSegment[z][y*width+x], considerIntensityCartesian.getValue(), component1, component2, component3, weight, distance1, width, height, depth, cpt);
						cpt++;
					}
				}


				// x
				// compute density for x distribution for each experiment
				float bandwidthForGaussianDistributionForX = compute_rule_of_thumb_bandwidth(component1);

				// compute x histograms
				float[] interpolatedXHistogram = new float[nbBinsForX.getValue()],
						xHistogram = computeSymmetricalHistogram(component1,weight,distance1,nbBinsForX.getValue(),minX,maxX,interpolatedXHistogram);

				// compute Gaussian kernel
				float[] GaussianKernelForX = computeSymmetricalGaussianKernel(nbBinsForX.getValue(),bandwidthForGaussianDistributionForX,maxX);

				// compute x density
				float[] xDensity = computeDensity(nbBinsForX.getValue(),xHistogram,GaussianKernelForX);

				// y
				// compute density for y distribution for each experiment
				float bandwidthForGaussianDistributionForY = compute_rule_of_thumb_bandwidth(component2);

				// compute x histograms
				float[] interpolatedYHistogram = new float[nbBinsForY.getValue()],
						yHistogram = computeSymmetricalHistogram(component2,weight,distance1,nbBinsForY.getValue(),minY,maxY,interpolatedYHistogram);

				// compute Gaussian kernel
				float[] GaussianKernelForY = computeSymmetricalGaussianKernel(nbBinsForY.getValue(),bandwidthForGaussianDistributionForY,maxY);

				// compute x density
				float[] yDensity = computeDensity(nbBinsForY.getValue(),yHistogram,GaussianKernelForY);

				// depth
				// compute density for depth distribution for each experiment
				if(depth>1){
					// compute density for depth distribution for each experiment
					float bandwidthForGaussianDistributionForDepth = compute_rule_of_thumb_bandwidth(component3);

					// compute depth histograms
					float[] interpolatedDepthHistogram = new float[nbBinsForCartesianDepth.getValue()],
							depthHistogram = computeDepthHistogram(component3,weight,distance2,nbBinsForCartesianDepth.getValue(),depth-1,interpolatedDepthHistogram);

					// compute Gaussian kernel for each experiment
					float[] GaussianKernelForDepth = computeGaussianKernel(nbBinsForCartesianDepth.getValue(),bandwidthForGaussianDistributionForDepth,depth-1);

					// compute depth density
					float[] depthDensity = computeDensity(nbBinsForCartesianDepth.getValue(),depthHistogram,GaussianKernelForDepth);

					File f = exportCartesianExcelFile.getValue(true);
					if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
					exportCartesianFiles(xHistogram,yHistogram,depthHistogram,xDensity,yDensity,depthDensity,f,inputImage);
				}
				else{
					File f = exportCartesianExcelFile.getValue(true);
					if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
					exportCartesianFiles(xHistogram,yHistogram,xDensity,yDensity,f,inputImage);
				}
			}

			else{
				// compute events coordinates
				if(coordinateSystem.getValue()=="Spherical"){

					if((ref1Xspherical.getValue()>-1)&&(ref1Xspherical.getValue()<width)&&(ref1Yspherical.getValue()>-1)&&(ref1Yspherical.getValue()<height)){
						referenceCenter1 = ref1Xspherical.getValue();
						referenceCenter2 = ref1Yspherical.getValue();
						if((ref1Zspherical.getValue()>-1)&&(ref1Zspherical.getValue()<depth)){
							referenceCenter3 = ref1Zspherical.getValue();
						}
						else{
							referenceCenter3 = cellCenter3;
						}
					}
					else{
						referenceCenter1 = cellCenter1;
						referenceCenter2 = cellCenter2;
						referenceCenter3 = cellCenter3;
					}
					// reference direction
					if((ref2Xspherical.getValue()>-1)&&(ref2Xspherical.getValue()<width)&&(ref2Yspherical.getValue()>-1)&&(ref2Yspherical.getValue()<height)){
						referenceDirection1 = (float)(Math.atan2((referenceCenter2-ref2Yspherical.getValue()),(ref2Xspherical.getValue()-referenceCenter1)) + Math.PI/2.);
					}
					else{
						if((referenceCenter1!=cellCenter1)&&(referenceCenter2!=cellCenter2)){
							referenceDirection1 = (float)(Math.atan2((referenceCenter2-cellCenter2),(cellCenter1-referenceCenter1)) + Math.PI/2.);
						}
					}

					float maxRadius=0;
					int nbTrajectories;
					component1 = new float[nbTrajectories];
					component2 = new float[nbTrajectories];
					component3 = new float[nbTrajectories];
					weight = new float[nbTrajectories];
					distance1 = new float[nbTrajectories];
					distance2 = new float[nbTrajectories];

					// compute normalizing distance
					// compute normalizing distances for each spatial point so it's done once for all
					float[][] cellSegment = new float[depth][arraySize];
					int[][] insideXCoord = new int[depth][arraySize];
					int[][] insideYCoord = new int[depth][arraySize];
					for(int z=0;z<depth;z++){
						for(int y=0;y<height;y++){
							for(int x=0;x<width;x++){
								if(intensitySumOverTime[y*width+x]>0.001){
									if((cellMaskImage!=null)||(forbiddenRegionImage!=null)){
										if(forbiddenRegionImage!=null){
											if(forbiddenRegionArray[z][referenceCenter2*width+referenceCenter1]>0){
												cellSegment[z][y*width+x] = computeDistanceToCellBorder(referenceCenter1,referenceCenter2,x,y,z,cellBorder,innerBorder,width,height,depth,insideXCoord,insideYCoord);
											}
											else{
												cellSegment[z][y*width+x] = computeSphericalDistanceToCellBorder(referenceCenter1,referenceCenter2,referenceCenter3,x,y,z,cellBorder,width,height,depth);
											}
										}
										else{
											cellSegment[z][y*width+x] = computeSphericalDistanceToCellBorder(referenceCenter1,referenceCenter2,referenceCenter3,x,y,z,cellBorder,width,height,depth);
										}
									}
									else{
										cellSegment[z][y*width+x] = 1;
									}

								}
							}
						}
					}

					int cpt=0;
					///////////////////////////////////////
					// loop over trajectories
					///////////////////////////////////////
					for(trajectories){
						// extract coordinates x,y,z associated to current trajectory
						// the middle point should be favoured
						if((cellMaskArray[z][y*width+x]>0)&&(inputArray[z][y*width+x]>0)){
							// compute angle
							float theta=(float)(Math.atan2(referenceCenter2-y,x-referenceCenter1));
							// take into account the orientation reference
							float orientation=(theta-referenceDirection1);
							while(orientation<0){orientation+=(2*Math.PI);}
							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

							// initialization
							float currentRadius=0;
							if(forbiddenRegionImage!=null){
								if(forbiddenRegionArray[z][referenceCenter2*width+referenceCenter1]>0){
									currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[z][y*width+x]-x,2.)+Math.pow(insideYCoord[z][y*width+x]-y,2.)));
								}
								else{
									currentRadius =	(float)(Math.sqrt(Math.pow(referenceCenter1-x,2.)+Math.pow(referenceCenter2-y,2.)));
								}
							}
							else{
								currentRadius =	(float)(Math.sqrt(Math.pow(referenceCenter1-x,2.)+Math.pow(referenceCenter2-y,2.)));
							}

							if((float)Math.sqrt(currentRadius*currentRadius+z*z)>maxRadius){maxRadius = (float)Math.sqrt(currentRadius*currentRadius+z*z);}
							extractSphericalCoordinates(inputArray, x, y, z, currentRadius, orientation, cellSegment[z][y*width+x], considerIntensitySpherical.getValue(), component1, component2, component3, weight, distance1, width, height, depth, cpt);
							cpt++;
						}
					}

					// radius
					// compute density for radius distribution for each experiment
					float bandwidthForGaussianDistributionForRadius = compute_rule_of_thumb_bandwidth(component1);

					// compute radius histograms
					// To modify distanceToZborder -> distanceToPlaneBorder
					float[] radiusHistogram = computeRadiusHistogram(component1,weight,distance1,nbBinsForSphericalRadius.getValue());

					// compute Gaussian kernel for each experiment
					float[] GaussianKernelForRadius = computeGaussianKernel(nbBinsForSphericalRadius.getValue(),bandwidthForGaussianDistributionForRadius,maxRadius);

					// compute radius density
					float[] radiusDensity = computeDensity(nbBinsForSphericalRadius.getValue(),radiusHistogram,GaussianKernelForRadius);


					// azimuthal angle
					// compute density for azimuthal angle distribution for each experiment
					float concentrationForVonMisesDistributionOfAzimuthalAngle = compute_rule_of_thumb_concentration(component2);

					// compute azimuthal angle histograms
					float[] azimuthalAngleHistogram = computeCircularHistogram(component2,weight,distance1,nbBinsForFirstAngle.getValue());

					// compute von Mises kernel for each experiment
					float[] vonMisesKernelOfAzimuthalAngle = computeVonMisesKernel(nbBinsForFirstAngle.getValue(),concentrationForVonMisesDistributionOfAzimuthalAngle);

					// compute azimuthal angle density
					float[] azimuthalAngleDensity = computeCircularDensity(nbBinsForFirstAngle.getValue(),azimuthalAngleHistogram,vonMisesKernelOfAzimuthalAngle);


					// polar angle
					// compute density for polar angle distribution for each experiment
					float concentrationForVonMisesDistributionOfPolarAngle = compute_rule_of_thumb_concentration(component3);

					// compute polar angle histograms
					float[] polarAngleHistogram = computeCircularHistogram(component3,weight,distance1,nbBinsForSecondAngle.getValue());

					// compute von Mises kernel for each experiment
					float[] vonMisesKernelOfPolarAngle = computeVonMisesKernel(nbBinsForSecondAngle.getValue(),concentrationForVonMisesDistributionOfPolarAngle);

					// compute polar angle density
					float[] polarAngleDensity = computeCircularDensity(nbBinsForSecondAngle.getValue(),polarAngleHistogram,vonMisesKernelOfPolarAngle);

					File f = exportSphericalExcelFile.getValue(true);
					if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
					exportSphericalFiles(radiusHistogram,azimuthalAngleHistogram,polarAngleHistogram,radiusDensity,azimuthalAngleDensity,polarAngleDensity,f,inputImage);
				}

				else{
					if(coordinateSystem.getValue()=="Distance to cell border"){

						float maxDistance=0;
						int nbTrajectories=0;
						for(int t=0;t<nbFrames;t++){
							float[][] inputArray = Array2DUtil.arrayToFloatArray(inputImage.getDataXYZ(t,channelOfInterestDensities.getValue()), inputImage.isSignedDataType());
							for(int z=0;z<depth;z++){
								for(int y=0;y<height;y++){
									for(int x=0;x<width;x++){
										if((cellMaskArray[z][y*width+x]>0)&&(forbiddenRegionArray[z][y*width+x]==0)&&(inputArray[z][y*width+x]>0)){
											nbTrajectories++;
										}
									}
								}
							}
						}
						component1 = new float[nbTrajectories];
						weight = new float[nbTrajectories];
						distance1 = new float[nbTrajectories];


						// compute normalizing distance
						float[][] cellSegment = new float[depth][arraySize];
						for(int z=0;z<depth;z++){
							for(int y=0;y<height;y++){
								for(int x=0;x<width;x++){
									if(intensitySumOverTime[y*width+x]>0.001){
										if((cellMaskImage!=null)||(forbiddenRegionImage!=null)){
											if(forbiddenRegionImage!=null){
												cellSegment[z][y*width+x] = computeActualDistanceToCellBorder(x,y,z,cellBorder,innerBorder,width,height,depth);
											}
											else{
												cellSegment[z][y*width+x] = computeActualDistanceToCellBorder(x,y,z,cellBorder,width,height,depth);
											}
										}
									}
								}
							}
						}

						int cpt=0;
						if(forbiddenRegionImage==null){
							for(trajectories){
								// extract coordinates x,y,z associated to current trajectory
								// the middle point should be favoured
								if((cellMaskArray[z][y*width+x]>0)&&(inputArray[z][y*width+x]>0)){
									if(cellSegment[z][y*width+x]>maxDistance){maxDistance = cellSegment[z][y*width+x];}
									extractDistanceCoordinates(inputArray, x, y, z, cellSegment[z][y*width+x], considerIntensityDistance.getValue(), component1, weight, distance1, width, height, depth, cpt);
									cpt++;
								}
							}
						}
						else{
							if((cellMaskArray[z][y*width+x]>0)&&(forbiddenRegionArray[z][y*width+x]==0)&&(inputArray[z][y*width+x]>0)){
								if(cellSegment[z][y*width+x]>maxDistance){maxDistance = cellSegment[z][y*width+x];}
								extractDistanceCoordinates(inputArray, x, y, z, cellSegment[z][y*width+x], considerIntensityDistance.getValue(), component1, weight, distance1, width, height, depth, cpt);
								cpt++;
							}
						}

						// distance to cell border
						// compute density for radius distribution for each experiment
						float bandwidthForGaussianDistributionForDistance = compute_rule_of_thumb_bandwidth(component1);

						// compute distance histograms
						float[] distanceHistogram = computeRadiusHistogram(component1,weight,distance1,nbBinsForDistance.getValue());

						// compute Gaussian kernel for each experiment
						float[] GaussianKernelForDistance = computeGaussianKernel(nbBinsForDistance.getValue(),bandwidthForGaussianDistributionForDistance,maxDistance);

						// compute distance density
						float[] distanceDensity = computeDensity(nbBinsForDistance.getValue(),distanceHistogram,GaussianKernelForDistance);

						File f = exportDistanceExcelFile.getValue(true);
						if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
						exportDistanceFiles(distanceHistogram,distanceDensity,f,inputImage);
					}
				}
			}

		}
	}

}