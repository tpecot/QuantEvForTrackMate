import java.io.File;

import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.concurrent.locks.ReentrantLock;

public class QuantEv{

	// parameters associated with cell mask
	int[][] computeParametersAssociatedWithCellMask(){
		int[][] cellMaskArray = new int[depth][arraySize];
		int nbPts=0;
		if(cellMaskImage.getSizeZ()==depth){
			cellMaskArray = Array2DUtil.arrayToIntArray(cellMaskImage.getDataXYZ(0,0), cellMaskImage.isSignedDataType());
			for(int z=0;z<cellMaskImage.getSizeZ();z++){
				for(int y=0;y<cellMaskImage.getSizeY();y++){
					for(int x=0;x<cellMaskImage.getSizeX();x++){
						if(cellMaskArray[z][y*width+x]>0){
							cellMaskArray[z][y*width+x] = 1;
							cellCenter1 += x;
							cellCenter2 += y;
							cellCenter3 += z;
							nbPts++;
						}
					}
				}
			}
			cellCenter1 = (int)((float)cellCenter1/(float)nbPts);
			cellCenter2 = (int)((float)cellCenter2/(float)nbPts);
			cellCenter3 = (int)((float)cellCenter3/(float)nbPts);
		}
		else{
			int[][] cellMaskArray2D = Array2DUtil.arrayToIntArray(cellMaskImage.getDataXYZ(0,0), cellMaskImage.isSignedDataType());
			for(int y=0;y<cellMaskImage.getSizeY();y++){
				for(int x=0;x<cellMaskImage.getSizeX();x++){
					if(cellMaskArray2D[0][y*width+x]>0){
						for(int z=0;z<depth;z++){
							cellMaskArray[z][y*width+x] = 1;
							cellCenter1 += x;
							cellCenter2 += y;
							nbPts++;
						}
					}
				}
			}
			cellCenter1 = (int)((float)cellCenter1/(float)nbPts);
			cellCenter2 = (int)((float)cellCenter2/(float)nbPts);
			cellCenter3 = depth/2;
		}
		return cellMaskArray;
	}

	// parameters associated with cell mask
	int[][] computeParametersAssociatedWithForbiddenRegionMask(){
		
		int[][] forbiddenRegionArray = new int[depth][arraySize];
		if(forbiddenRegionImage.getSizeZ()==depth){
			forbiddenRegionArray = Array2DUtil.arrayToIntArray(forbiddenRegionImage.getDataXYZ(0,0), forbiddenRegionImage.isSignedDataType());
			for(int z=0;z<depth;z++){
				for(int y=0;y<height;y++){
					for(int x=0;x<width;x++){
						if(forbiddenRegionArray[z][y*width+x]>0){
							forbiddenRegionArray[z][y*width+x] = 1;
						}
					}
				}
			}
		}
		else{
			int [][] forbiddenRegionArray2D = Array2DUtil.arrayToIntArray(forbiddenRegionImage.getDataXYZ(0,0), forbiddenRegionImage.isSignedDataType());
			for(int y=0;y<height;y++){
				for(int x=0;x<width;x++){
					if(forbiddenRegionArray2D[0][y*width+x]>0){
						for(int z=0;z<depth;z++){
							forbiddenRegionArray[z][y*width+x] = 1;
						}
					}
				}
			}
		}
		
		return forbiddenRegionArray;
	}
	
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
// parameters associated with cell mask
	int[][] computeParametersAssociatedWithCellMask(){
		int[][] cellMaskArray = new int[depth][arraySize];
		int nbPts=0;
		if(cellMaskImage.getSizeZ()==depth){
			cellMaskArray = Array2DUtil.arrayToIntArray(cellMaskImage.getDataXYZ(0,0), cellMaskImage.isSignedDataType());
			for(int z=0;z<cellMaskImage.getSizeZ();z++){
				for(int y=0;y<cellMaskImage.getSizeY();y++){
					for(int x=0;x<cellMaskImage.getSizeX();x++){
						if(cellMaskArray[z][y*width+x]>0){
							cellMaskArray[z][y*width+x] = 1;
							cellCenter1 += x;
							cellCenter2 += y;
							cellCenter3 += z;
							nbPts++;
						}
					}
				}
			}
			cellCenter1 = (int)((float)cellCenter1/(float)nbPts);
			cellCenter2 = (int)((float)cellCenter2/(float)nbPts);
			cellCenter3 = (int)((float)cellCenter3/(float)nbPts);
		}
		else{
			int[][] cellMaskArray2D = Array2DUtil.arrayToIntArray(cellMaskImage.getDataXYZ(0,0), cellMaskImage.isSignedDataType());
			for(int y=0;y<cellMaskImage.getSizeY();y++){
				for(int x=0;x<cellMaskImage.getSizeX();x++){
					if(cellMaskArray2D[0][y*width+x]>0){
						for(int z=0;z<depth;z++){
							cellMaskArray[z][y*width+x] = 1;
							cellCenter1 += x;
							cellCenter2 += y;
							nbPts++;
						}
					}
				}
			}
			cellCenter1 = (int)((float)cellCenter1/(float)nbPts);
			cellCenter2 = (int)((float)cellCenter2/(float)nbPts);
			cellCenter3 = depth/2;
		}
		return cellMaskArray;
	}

	// parameters associated with cell mask
	int[][] computeParametersAssociatedWithForbiddenRegionMask(){
		
		int[][] forbiddenRegionArray = new int[depth][arraySize];
		if(forbiddenRegionImage.getSizeZ()==depth){
			forbiddenRegionArray = Array2DUtil.arrayToIntArray(forbiddenRegionImage.getDataXYZ(0,0), forbiddenRegionImage.isSignedDataType());
			for(int z=0;z<depth;z++){
				for(int y=0;y<height;y++){
					for(int x=0;x<width;x++){
						if(forbiddenRegionArray[z][y*width+x]>0){
							forbiddenRegionArray[z][y*width+x] = 1;
						}
					}
				}
			}
		}
		else{
			int [][] forbiddenRegionArray2D = Array2DUtil.arrayToIntArray(forbiddenRegionImage.getDataXYZ(0,0), forbiddenRegionImage.isSignedDataType());
			for(int y=0;y<height;y++){
				for(int x=0;x<width;x++){
					if(forbiddenRegionArray2D[0][y*width+x]>0){
						for(int z=0;z<depth;z++){
							forbiddenRegionArray[z][y*width+x] = 1;
						}
					}
				}
			}
		}
		
		return forbiddenRegionArray;
	}

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
	void extractCylindricalCoordinates(int x,int y,int z,float rad,float orientation,float segment,float trackFeature,
			float[] radius,float[] angle,float[] depth,float[] intensity,float[] distanceToBorder,
			int dimX,int dimY,int dimZ,int i){
		angle[i] = orientation;
		intensity[i] = trackFeature;
		depth[i] = z;
		radius[i] = rad;
		distanceToBorder[i] = segment;
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
	void extractSphericalCoordinates(int x,int y,int z,float rad,float orientation,float segment,float trackFeature,
			float[] radius,float[] angle1,float[] angle2,float[] intensity,float[] distanceToBorder,
			int dimX,int dimY,int dimZ,int i){
		angle1[i] = orientation;
		intensity[i] = trackFeature;
		angle2[i] = (float)Math.acos(z/Math.sqrt(rad*rad+z*z));
		radius[i] = (float)Math.sqrt(rad*rad+z*z);
		distanceToBorder[i] = segment;
	}

	// extract event coordinates from images
	void extractDistanceCoordinates(int x,int y,int z,float segment,float trackFeature,
			float[] distance,float[] intensity,float[] unity,int dimX,int dimY,int dimZ,int i){
		distance[i] = segment;
		intensity[i] = trackFeature;
		unity[i] = 1;
	}


	// extract event coordinates from images
	void extractCartesianCoordinates(int x,int y,int z,float rad,float orientation,float segment,float trackFeature,
			float[] xTab,float[] yTab,float[] depth,float[] intensity,float[] distanceToBorder,
			int dimX,int dimY,int dimZ,int i){

		xTab[i] = rad*(float)Math.cos(orientation);
		yTab[i] = rad*(float)Math.sin(orientation);
		intensity[i] = trackFeature;
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
	float[] computeSymmetricalHistogram(float[] input,float[] weight,float[] length,int nbBins,float[] interpolatedHistogram){
		float minValue=100000,maxValue=0;
		for(int i=0;i<input.length;i++){
			if(!Float.isNaN(input[i])){
				if(input[i]<minValue){minValue = input[i];}
				if(input[i]>maxValue){maxValue = input[i];}
			}
		}
		
		int actualNbBins;
		if(nbBins>(2*maxValue+1.)){actualNbBins = (int)Math.ceil(2*maxValue+1);}
		else{actualNbBins = nbBins;}
		float[] histogram = new float[actualNbBins],
				outputHistogram = new float[nbBins];
		int[] nbEltsPerBin = new int[nbBins];
		for(int i=0;i<input.length;i++){
			if(!Float.isNaN(input[i])){
				int currentBin = Math.round((float)(actualNbBins-1)*(input[i]-minValue)/(2*maxValue));
				if(currentBin>=actualNbBins){currentBin=actualNbBins-1;}
				if(currentBin<0){currentBin=0;}
				if(!Float.isNaN(weight[i])&&!Float.isNaN(length[i])&&(weight[i]<1000000)){
					histogram[currentBin] += weight[i];
					nbEltsPerBin[currentBin]++;
				}
			}
		}
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
	float[] computeSymmetricalGaussianKernel(int nbBins,float bandwidth,float[] input){
		float maxValue=0;
		for(int i=0;i<input.length;i++){
			if(!Float.isNaN(input[i])){
				if(input[i]>maxValue){maxValue = input[i];}
			}
		}
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
	//return the concentration for the van mises distribution
	//param S a vector of angles
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
		for(int i=0;i<nbBins;i++){
			float vonMisesNormalization=0;
			for(int j=0;j<nbBins;j++){
				int currentIndex = i+j-nbBins/2;
				if(currentIndex<0){currentIndex += nbBins;}
				if(currentIndex>=nbBins){currentIndex -= nbBins;}
				density[i] += histogram[currentIndex]*vonMisesKernel[j];
				vonMisesNormalization += vonMisesKernel[j];
			}
			density[i] = density[i]/vonMisesNormalization;
		}

		return density;
	}

	File getValidSaveFile(){
        boolean hasValidSaveFile = false;
        File file = null;
        while (!hasValidSaveFile)
        {
            JFileChooser fileChooser = new JFileChooser();
            fileChooser.setMultiSelectionEnabled(false);
            fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
            fileChooser.setName("Text file for saving results");
            int returnVal = fileChooser.showDialog(panel, "Set as save file");
            if (returnVal == JFileChooser.APPROVE_OPTION)
            {
                file = fileChooser.getSelectedFile();
                if (file.exists())
                {
                    int n = JOptionPane.showConfirmDialog(Icy.getMainInterface().getMainFrame(),
                            "This file already exists. Do you want to add QuantEv analysis in this same file?",
                            "Save QuantEv results", JOptionPane.YES_NO_CANCEL_OPTION);
                    switch (n)
                    {
                        case JOptionPane.YES_OPTION:
                            hasValidSaveFile = true;
                            break;
                        case JOptionPane.CANCEL_OPTION:
                            hasValidSaveFile = true;
                            file = null;
                            break;
                        case JOptionPane.NO_OPTION:
                            hasValidSaveFile = false;
                    }
                }
                else
                    hasValidSaveFile = true;
            }
            else
                return null;
        }
        return file;
    }

	// export histograms and densities as xls file
	void exportCylindricalFiles(float[] histogram1,float[] histogram2,float[] histogram3,float[] density1,float[] density2,float[] density3,File f){
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

			if(histogram1.length!=(radiusRows)){MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, currentColumn, bin, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, currentColumn, bin, density1[bin]);
			}
			if(histogram2.length!=(angleRows)){MessageDialog.showDialog("In the excel file, the worksheets Polar angle histogram and Polar angle density do not have the same number of rows than the number of bins for polar angle.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, 0, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, currentColumn, 0, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, bin, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, currentColumn, bin, density2[histogram2.length-bin]);
			}
			if(histogram3.length!=(depthRows)){MessageDialog.showDialog("In the excel file, the worksheets Depth histogram and Depth density do not have the same number of rows than the number of bins for depth.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(depthHistogramSheet, currentColumn, bin, histogram3[bin]);
				XLSUtil.setCellNumber(depthDensitySheet, currentColumn, bin, density3[bin]);
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
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, 0, bin, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, 0, bin, density1[bin]);
			}
			XLSUtil.setCellNumber(angleHistogramSheet, 0, 0, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, 0, 0, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, 0, bin, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, 0, bin, density2[histogram2.length-bin]);
			}
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(depthHistogramSheet, 0, bin, histogram3[bin]);
				XLSUtil.setCellNumber(depthDensitySheet, 0, bin, density3[bin]);
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
	void exportCylindricalFiles(float[] histogram1,float[] histogram2,float[] density1,float[] density2,File f){
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

			if(histogram1.length!=(radiusRows)){MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, currentColumn, bin, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, currentColumn, bin, density1[bin]);
			}
			if(histogram2.length!=(angleRows)){MessageDialog.showDialog("In the excel file, the worksheets Polar angle histogram and Polar angle density do not have the same number of rows than the number of bins for polar angle.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, 0, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, currentColumn, 0, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, currentColumn, bin, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, currentColumn, bin, density2[histogram2.length-bin]);
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
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, 0, bin, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, 0, bin, density1[bin]);
			}
			XLSUtil.setCellNumber(angleHistogramSheet, 0, 0, histogram2[0]);
			XLSUtil.setCellNumber(angleDensitySheet, 0, 0, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angleHistogramSheet, 0, bin, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angleDensitySheet, 0, bin, density2[histogram2.length-bin]);
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
	void exportCartesianFiles(float[] histogram1,float[] histogram2,float[] histogram3,float[] density1,float[] density2,float[] density3,File f){
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

			if(histogram1.length!=(xRows)){MessageDialog.showDialog("In the excel file, the worksheets X histogram and X density do not have the same number of rows than the number of bins for abscissa.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, currentColumn, bin, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, currentColumn, bin, density1[bin]);
			}
			if(histogram2.length!=(yRows)){MessageDialog.showDialog("In the excel file, the worksheets Y histogram and Y density do not have the same number of rows than the number of bins for ordinate.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, currentColumn, bin, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, currentColumn, bin, density2[bin]);
			}
			if(histogram3.length!=(zRows)){MessageDialog.showDialog("In the excel file, the worksheets Z histogram and Z density do not have the same number of rows than the number of bins for height.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(zHistogramSheet, currentColumn, bin, histogram3[bin]);
				XLSUtil.setCellNumber(zDensitySheet, currentColumn, bin, density3[bin]);
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
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, 0, bin, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, 0, bin, density1[bin]);
			}
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, 0, bin, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, 0, bin, density2[bin]);
			}
			for(int bin=0;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(zHistogramSheet, 0, bin, histogram3[bin]);
				XLSUtil.setCellNumber(zDensitySheet, 0, bin, density3[bin]);
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
	void exportCartesianFiles(float[] histogram1,float[] histogram2,float[] density1,float[] density2,File f){
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

			if(histogram1.length!=(xRows)){MessageDialog.showDialog("In the excel file, the worksheets X histogram and X density do not have the same number of rows than the number of bins for abscissa.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, currentColumn, bin, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, currentColumn, bin, density1[bin]);
			}
			if(histogram2.length!=(yRows)){MessageDialog.showDialog("In the excel file, the worksheets Y histogram and Y density do not have the same number of rows than the number of bins for ordinate.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, currentColumn, bin, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, currentColumn, bin, density2[bin]);
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
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(xHistogramSheet, 0, bin, histogram1[bin]);
				XLSUtil.setCellNumber(xDensitySheet, 0, bin, density1[bin]);
			}
			for(int bin=0;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(yHistogramSheet, 0, bin, histogram2[bin]);
				XLSUtil.setCellNumber(yDensitySheet, 0, bin, density2[bin]);
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
	void exportSphericalFiles(float[] histogram1,float[] histogram2,float[] histogram3,float[] density1,float[] density2,float[] density3,File f){
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

			if(histogram1.length!=(radiusRows)){MessageDialog.showDialog("In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, currentColumn, bin, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, currentColumn, bin, density1[bin]);
			}
			if(histogram2.length!=(angle1Rows)){MessageDialog.showDialog("In the excel file, the worksheets Colatitude histogram and Colatitude density do not have the same number of rows than the number of bins for colatitude.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellNumber(angle1HistogramSheet, currentColumn, 0, histogram2[0]);
			XLSUtil.setCellNumber(angle1DensitySheet, currentColumn, 0, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angle1HistogramSheet, currentColumn, bin, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angle1DensitySheet, currentColumn, bin, density2[histogram2.length-bin]);
			}
			if(histogram3.length!=(angle2Rows)){MessageDialog.showDialog("In the excel file, the worksheets Azimuth angle histogram and Azimuth angle density do not have the same number of rows than the number of bins for azimuth angle.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			XLSUtil.setCellNumber(angle2HistogramSheet, currentColumn, 0, histogram3[0]);
			XLSUtil.setCellNumber(angle2DensitySheet, currentColumn, 0, density3[0]);
			for(int bin=1;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(angle2HistogramSheet, currentColumn, bin, histogram2[histogram3.length-bin]);
				XLSUtil.setCellNumber(angle2DensitySheet, currentColumn, bin, density2[histogram3.length-bin]);
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
			for(int bin=0;bin<histogram1.length;bin++){
				XLSUtil.setCellNumber(radiusHistogramSheet, 0, bin, histogram1[bin]);
				XLSUtil.setCellNumber(radiusDensitySheet, 0, bin, density1[bin]);
			}
			XLSUtil.setCellNumber(angle1HistogramSheet, 0, 0, histogram2[0]);
			XLSUtil.setCellNumber(angle1DensitySheet, 0, 0, density2[0]);
			for(int bin=1;bin<histogram2.length;bin++){
				XLSUtil.setCellNumber(angle1HistogramSheet, 0, bin, histogram2[histogram2.length-bin]);
				XLSUtil.setCellNumber(angle1DensitySheet, 0, bin, density2[histogram2.length-bin]);
			}
			XLSUtil.setCellNumber(angle2HistogramSheet, 0, 0, histogram3[0]);
			XLSUtil.setCellNumber(angle2DensitySheet, 0, 0, density3[0]);
			for(int bin=1;bin<histogram3.length;bin++){
				XLSUtil.setCellNumber(angle2HistogramSheet, 0, bin, histogram3[histogram3.length-bin]);
				XLSUtil.setCellNumber(angle2DensitySheet, 0, bin, density3[histogram3.length-bin]);
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
	void exportDistanceFiles(float[] histogram,float[] density,File f){
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

			if(histogram.length!=(distanceRows)){MessageDialog.showDialog("In the excel file, the worksheets Distance to cell border histogram and Distance to cell border density do not have the same number of rows than the number of bins for distance to cell border.");
			try
			{
				XLSUtil.saveAndClose(wb);
			}
			catch (Exception e)
			{
				throw new IcyHandledException(e.getMessage());
			}
			return;}
			for(int bin=0;bin<histogram.length;bin++){
				XLSUtil.setCellNumber(distanceHistogramSheet, currentColumn, bin, histogram[bin]);
				XLSUtil.setCellNumber(distanceDensitySheet, currentColumn, bin, density[bin]);
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
			for(int bin=0;bin<histogram.length;bin++){
				XLSUtil.setCellNumber(distanceHistogramSheet, 0, bin, histogram[bin]);
				XLSUtil.setCellNumber(distanceDensitySheet, 0, bin, density[bin]);
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

	// compute track feature
	float computeTrackFeature(TrackSegment ts,String feat){
		float feature=0;
		if(feat=="Confinement ratio"){
			Detection dFirst = ts.getDetectionList().get(0), dLast = ts.getDetectionList().get(ts.getDetectionList().size()-1);
			float dd = (float)Math.sqrt(Math.pow(dFirst.getX()-dLast.getX(),2.)+Math.pow(dFirst.getY()-dLast.getY(),2.)+Math.pow(dFirst.getZ()-dLast.getZ(),2.));
			float tpl=0;
			for(int i=0;i<(ts.getDetectionList().size()-1);i++){
				Detection currentPt = ts.getDetectionList().get(i),nextPt = ts.getDetectionList().get(i+1);
				tpl += (float)Math.sqrt(Math.pow(nextPt.getX()-currentPt.getX(),2.)+Math.pow(nextPt.getY()-currentPt.getY(),2.)+Math.pow(nextPt.getZ()-currentPt.getZ(),2.)); 
			}
			if(tpl>0.001){
				feature = dd/tpl;
			}
		}
		if(feat=="Displacement distance"){
			Detection dFirst = ts.getDetectionList().get(0), dLast = ts.getDetectionList().get(ts.getDetectionList().size()-1);
			feature = (float)Math.sqrt(Math.pow(dFirst.getX()-dLast.getX(),2.)+Math.pow(dFirst.getY()-dLast.getY(),2.)+Math.pow(dFirst.getZ()-dLast.getZ(),2.));
		}
		if(feat=="Total path length"){
			for(int i=0;i<(ts.getDetectionList().size()-1);i++){
				Detection currentPt = ts.getDetectionList().get(i),nextPt = ts.getDetectionList().get(i+1);
				feature += (float)Math.sqrt(Math.pow(nextPt.getX()-currentPt.getX(),2.)+Math.pow(nextPt.getY()-currentPt.getY(),2.)+Math.pow(nextPt.getZ()-currentPt.getZ(),2.)); 
			}
		}
		if(feat=="Lifetime"){
			feature = (float)ts.getDetectionList().size();
		}
		
		return feature;
	}

	JComboBox<String> InputTypeBox,TrackFeatureBox;
	JPanel CylindricalPanel,SphericalPanel,CartesianPanel,DistancePanel;
	JPanel CoordinateSystemPanel;
	CardLayout CoordinateSystemLayout;

	SpinnerNumberModel nbBinsForX = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForY = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForCartesianDepth = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForCylindricalRadius = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForAngle = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForCylindricalDepth = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForSphericalRadius = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForFirstAngle = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForSecondAngle = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForDistance = new SpinnerNumberModel(30,5,1000,1),
			ref1Xcartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Ycartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Xcylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Ycylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Xspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Yspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Zspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Xcartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Ycartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Xcylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Ycylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Xspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Yspherical = new SpinnerNumberModel(-1,-1,100000,1);
	
	Sequence cellMaskImage = new Sequence("Cell mask"),forbiddenRegionImage = new Sequence("Exclusion mask");
	int[][] cellMaskArray,forbiddenRegionArray;
	
	int arraySize=0,width=0,height=0,depth=0,nbFrames=0;
	int cellCenter1=0, cellCenter2=0, cellCenter3=0;
	int minX=10000,minY=10000,minZ=10000,maxX=0,maxY=0,maxZ=0;
	
				
				
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
	
	JComboBox<String> InputTypeBox,TrackFeatureBox;
	JPanel CylindricalPanel,SphericalPanel,CartesianPanel,DistancePanel;
	JPanel CoordinateSystemPanel;
	CardLayout CoordinateSystemLayout;

	SpinnerNumberModel nbBinsForX = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForY = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForCartesianDepth = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForCylindricalRadius = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForAngle = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForCylindricalDepth = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForSphericalRadius = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForFirstAngle = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForSecondAngle = new SpinnerNumberModel(30,5,1000,1),
			nbBinsForDistance = new SpinnerNumberModel(30,5,1000,1),
			ref1Xcartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Ycartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Xcylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Ycylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Xspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Yspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref1Zspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Xcartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Ycartesian = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Xcylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Ycylindrical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Xspherical = new SpinnerNumberModel(-1,-1,100000,1),
			ref2Yspherical = new SpinnerNumberModel(-1,-1,100000,1);
	
	Sequence cellMaskImage = new Sequence("Cell mask"),forbiddenRegionImage = new Sequence("Exclusion mask");
	int[][] cellMaskArray,forbiddenRegionArray;
	
	int arraySize=0,width=0,height=0,depth=0,nbFrames=0;
	int cellCenter1=0, cellCenter2=0, cellCenter3=0;
	int minX=10000,minY=10000,minZ=10000,maxX=0,maxY=0,maxZ=0;
			
	public QuantEvTrackProcessor()
	{
		this.setName("QuantEv");
		
	}

	@Override
	public void Close()
	{

	}

	@Override
	public void Compute()
	{
		this.panel.setLayout(new BorderLayout());
		
		JPanel northPanel = new JPanel();
		northPanel.setLayout(new BorderLayout());

		CoordinateSystemLayout = new CardLayout();
		CoordinateSystemPanel = new JPanel(CoordinateSystemLayout);	
		
		CylindricalPanel = new JPanel(new GridLayout(9, 2));
		CylindricalPanel.add(new JLabel("X coordinate of cylindrical coordinate system center"));
		JSpinner ref1XcylindricalSpinner = new JSpinner(ref1Xcylindrical);
		CylindricalPanel.add(ref1XcylindricalSpinner);
		CylindricalPanel.add(new JLabel("Y coordinate of cylindrical coordinate system center"));
		JSpinner ref1YcylindricalSpinner = new JSpinner(ref1Ycylindrical);
		CylindricalPanel.add(ref1YcylindricalSpinner);
		CylindricalPanel.add(new JLabel("X coordinate of point located on the polar axis"));
		JSpinner ref2XcylindricalSpinner = new JSpinner(ref2Xcylindrical);
		CylindricalPanel.add(ref2XcylindricalSpinner);
		CylindricalPanel.add(new JLabel("Y coordinate of point located on the polar axis"));
		JSpinner ref2YcylindricalSpinner = new JSpinner(ref2Ycylindrical);
		CylindricalPanel.add(ref2YcylindricalSpinner);
		CylindricalPanel.add(new JLabel("Number of bins for radius"));
		JSpinner nbBinsForCylindricalRadiusSpinner = new JSpinner(nbBinsForCylindricalRadius);
		CylindricalPanel.add(nbBinsForCylindricalRadiusSpinner);
		CylindricalPanel.add(new JLabel("Number of bins for polar angle"));
		JSpinner nbBinsForAngleSpinner = new JSpinner(nbBinsForAngle);
		CylindricalPanel.add(nbBinsForAngleSpinner);
		CylindricalPanel.add(new JLabel("Number of bins for depth"));
		JSpinner nbBinsForCylindricalDepthSpinner = new JSpinner(nbBinsForCylindricalDepth);
		CylindricalPanel.add(nbBinsForCylindricalDepthSpinner);
		CylindricalPanel.add(new JLabel("Cell mask"));
		SequenceChooser cellMaskCylindrical = new SequenceChooser();
		cellMaskCylindrical.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					cellMaskImage = cellMaskCylindrical.getSelectedSequence();
					arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
					width=cellMaskImage.getSizeX();
					height=cellMaskImage.getSizeY();
					depth=maxZ+1;
					cellMaskArray = computeParametersAssociatedWithCellMask();
				} catch (Exception eCMC) {}
			}
		});
		CylindricalPanel.add(cellMaskCylindrical);
		CylindricalPanel.add(new JLabel("Exclusion mask"));
		SequenceChooser forbiddenRegionCylindrical = new SequenceChooser();
		forbiddenRegionCylindrical.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					forbiddenRegionImage = forbiddenRegionCylindrical.getSelectedSequence();
					arraySize=forbiddenRegionImage.getSizeX()*forbiddenRegionImage.getSizeY();
					width=forbiddenRegionImage.getSizeX();
					height=forbiddenRegionImage.getSizeY();
					depth=maxZ+1;
					forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
				} catch (Exception eFRC) {}
			}
		});
		CylindricalPanel.add(forbiddenRegionCylindrical);
		CoordinateSystemPanel.add(CylindricalPanel,"Cylindrical");
		
		SphericalPanel = new JPanel(new GridLayout(10, 2));
		SphericalPanel.add(new JLabel("X coordinate of spherical coordinate system center"));
		JSpinner ref1XsphericalSpinner = new JSpinner(ref1Xspherical);
		SphericalPanel.add(ref1XsphericalSpinner);
		SphericalPanel.add(new JLabel("Y coordinate of spherical coordinate system center"));
		JSpinner ref1YsphericalSpinner = new JSpinner(ref1Yspherical);
		SphericalPanel.add(ref1YsphericalSpinner);
		SphericalPanel.add(new JLabel("Z coordinate of spherical coordinate system center"));
		JSpinner ref1ZsphericalSpinner = new JSpinner(ref1Zspherical);
		SphericalPanel.add(ref1ZsphericalSpinner);
		SphericalPanel.add(new JLabel("X coordinate of point located on the colatitude axis"));
		JSpinner ref2XsphericalSpinner = new JSpinner(ref2Xspherical);
		SphericalPanel.add(ref2XsphericalSpinner);
		SphericalPanel.add(new JLabel("Y coordinate of point located on the colatitude axis"));
		JSpinner ref2YsphericalSpinner = new JSpinner(ref2Yspherical);
		SphericalPanel.add(ref2YsphericalSpinner);
		SphericalPanel.add(new JLabel("Number of bins for radial distance"));
		JSpinner nbBinsForSphericalRadiusSpinner = new JSpinner(nbBinsForSphericalRadius);
		SphericalPanel.add(nbBinsForSphericalRadiusSpinner);
		SphericalPanel.add(new JLabel("Number of bins for colatitude"));
		JSpinner nbBinsForFirstAngleSpinner = new JSpinner(nbBinsForFirstAngle);
		SphericalPanel.add(nbBinsForFirstAngleSpinner);
		SphericalPanel.add(new JLabel("Number of bins for azimuth angle"));
		JSpinner nbBinsForSecondAngleSpinner = new JSpinner(nbBinsForSecondAngle);
		SphericalPanel.add(nbBinsForSecondAngleSpinner);
		SphericalPanel.add(new JLabel("Cell mask"));
		SequenceChooser cellMaskSpherical = new SequenceChooser();
		cellMaskSpherical.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					cellMaskImage = cellMaskSpherical.getSelectedSequence();
					arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
					width=cellMaskImage.getSizeX();
					height=cellMaskImage.getSizeY();
					depth=maxZ+1;
					cellMaskArray = computeParametersAssociatedWithCellMask();
				} catch (Exception eCMS) {}
			}
		});
		SphericalPanel.add(cellMaskSpherical);
		SphericalPanel.add(new JLabel("Exclusion mask"));
		SequenceChooser forbiddenRegionSpherical = new SequenceChooser();
		forbiddenRegionSpherical.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					forbiddenRegionImage = forbiddenRegionSpherical.getSelectedSequence();
					arraySize=forbiddenRegionImage.getSizeX()*forbiddenRegionImage.getSizeY();
					width=forbiddenRegionImage.getSizeX();
					height=forbiddenRegionImage.getSizeY();
					depth=maxZ+1;
					forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
				} catch (Exception eFRS) {}
			}
		});
		SphericalPanel.add(forbiddenRegionSpherical);
		CoordinateSystemPanel.add(SphericalPanel,"Spherical");
		
		CartesianPanel = new JPanel(new GridLayout(9, 2));
		CartesianPanel.add(new JLabel("X coordinate of Cartesian coordinate system center"));
		JSpinner ref1XcartesianSpinner = new JSpinner(ref1Xcartesian);
		CartesianPanel.add(ref1XcartesianSpinner);
		CartesianPanel.add(new JLabel("X coordinate of Cartesian coordinate system center"));
		JSpinner ref1YcartesianSpinner = new JSpinner(ref1Ycartesian);
		CartesianPanel.add(ref1YcartesianSpinner);
		CartesianPanel.add(new JLabel("X coordinate of point located on the X axis"));
		JSpinner ref2XcartesianSpinner = new JSpinner(ref2Xcartesian);
		CartesianPanel.add(ref2XcartesianSpinner);
		CartesianPanel.add(new JLabel("Y coordinate of point located on the X axis"));
		JSpinner ref2YcartesianSpinner = new JSpinner(ref2Ycartesian);
		CartesianPanel.add(ref2YcartesianSpinner);
		CartesianPanel.add(new JLabel("Number of bins for abscissa"));
		JSpinner nbBinsForXSpinner = new JSpinner(nbBinsForX);
		CartesianPanel.add(nbBinsForXSpinner);
		CartesianPanel.add(new JLabel("Number of bins for ordinate"));
		JSpinner nbBinsForYSpinner = new JSpinner(nbBinsForY);
		CartesianPanel.add(nbBinsForYSpinner);
		CartesianPanel.add(new JLabel("Number of bins for height"));
		JSpinner nbBinsForCartesianDepthSpinner = new JSpinner(nbBinsForCartesianDepth);
		CartesianPanel.add(nbBinsForCartesianDepthSpinner);
		CartesianPanel.add(new JLabel("Cell mask"));
		SequenceChooser cellMaskCartesian = new SequenceChooser();
		cellMaskCartesian.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					cellMaskImage = cellMaskCartesian.getSelectedSequence();
					arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
					width=cellMaskImage.getSizeX();
					height=cellMaskImage.getSizeY();
					depth=maxZ+1;
					cellMaskArray = computeParametersAssociatedWithCellMask();
				} catch (Exception eCMCa) {}
			}
		});
		CartesianPanel.add(cellMaskCartesian);
		CartesianPanel.add(new JLabel("Exclusion mask"));
		SequenceChooser forbiddenRegionCartesian = new SequenceChooser();
		forbiddenRegionCartesian.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					forbiddenRegionImage = forbiddenRegionCartesian.getSelectedSequence();
					arraySize=forbiddenRegionImage.getSizeX()*forbiddenRegionImage.getSizeY();
					width=forbiddenRegionImage.getSizeX();
					height=forbiddenRegionImage.getSizeY();
					depth=maxZ+1;
					forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
				} catch (Exception eFRCa) {}
			}
		});
		CartesianPanel.add(forbiddenRegionCartesian);
		CoordinateSystemPanel.add(CartesianPanel,"Cartesian");
		
		DistancePanel = new JPanel(new GridLayout(3, 2));
		DistancePanel.add(new JLabel("Number of bins for distance to cell border"));
		JSpinner nbBinsForDistanceSpinner = new JSpinner(nbBinsForDistance);
		DistancePanel.add(nbBinsForDistanceSpinner);
		DistancePanel.add(new JLabel("Cell mask"));
		SequenceChooser cellMaskDistance = new SequenceChooser();
		cellMaskDistance.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					cellMaskImage = cellMaskDistance.getSelectedSequence();
					arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
					width=cellMaskImage.getSizeX();
					height=cellMaskImage.getSizeY();
					depth=maxZ+1;
					cellMaskArray = computeParametersAssociatedWithCellMask();
				} catch (Exception eCMD) {}
			}
		});
		DistancePanel.add(cellMaskDistance);
		DistancePanel.add(new JLabel("Exclusion mask"));
		SequenceChooser forbiddenRegionDistance = new SequenceChooser();
		forbiddenRegionDistance.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				try{
					forbiddenRegionImage = forbiddenRegionDistance.getSelectedSequence();
					arraySize=forbiddenRegionImage.getSizeX()*forbiddenRegionImage.getSizeY();
					width=forbiddenRegionImage.getSizeX();
					height=forbiddenRegionImage.getSizeY();
					depth=maxZ+1;
					forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
				} catch (Exception eFRD) {}
			}
		});
		DistancePanel.add(forbiddenRegionDistance);
		CoordinateSystemPanel.add(DistancePanel,"Distance");
		
		String[] coordinateSystemPossibilities = {"Cylindrical","Spherical","Cartesian","Distance to cell border"};
		JPanel InputChoicePanel = new JPanel(new GridLayout(4, 1));
		InputChoicePanel.add(new JLabel("Coordinate system choice"));
		InputTypeBox = new JComboBox<String>(coordinateSystemPossibilities);
		InputChoicePanel.add(InputTypeBox);
		InputTypeBox.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(InputTypeBox.getSelectedItem().toString()=="Cylindrical"){
					CoordinateSystemLayout.show(CoordinateSystemPanel, "Cylindrical");
					try{
						cellMaskCylindrical.getSelectedItem().toString();
						cellMaskImage = cellMaskSpherical.getSelectedSequence();
						arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
						width=cellMaskImage.getSizeX();
						height=cellMaskImage.getSizeY();
						depth=maxZ+1;
						cellMaskArray = computeParametersAssociatedWithCellMask();
					} catch (Exception eFRD) {
						cellMaskArray = new int[0][0];
					}
					try{
						forbiddenRegionImage = forbiddenRegionCylindrical.getSelectedSequence();
						forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
					} catch (Exception eFRD) {
						forbiddenRegionArray = new int[0][0];
					}
				}
				if(InputTypeBox.getSelectedItem().toString()=="Spherical"){
					CoordinateSystemLayout.show(CoordinateSystemPanel, "Spherical");
					try{
						cellMaskSpherical.getSelectedItem().toString();
						cellMaskImage = cellMaskSpherical.getSelectedSequence();
						arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
						width=cellMaskImage.getSizeX();
						height=cellMaskImage.getSizeY();
						depth=maxZ+1;
						cellMaskArray = computeParametersAssociatedWithCellMask();
					} catch (Exception eFRD) {
						cellMaskArray = new int[0][0];
					}
					try{
						forbiddenRegionImage = forbiddenRegionSpherical.getSelectedSequence();
						forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
					} catch (Exception eFRD) {
						forbiddenRegionArray = new int[0][0];
					}
				}
				if(InputTypeBox.getSelectedItem().toString()=="Cartesian"){
					CoordinateSystemLayout.show(CoordinateSystemPanel, "Cartesian");
					try{
						cellMaskCartesian.getSelectedItem().toString();
						cellMaskImage = cellMaskCartesian.getSelectedSequence();
						arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
						width=cellMaskImage.getSizeX();
						height=cellMaskImage.getSizeY();
						depth=maxZ+1;
						cellMaskArray = computeParametersAssociatedWithCellMask();
					} catch (Exception eFRD) {
						cellMaskArray = new int[0][0];
					}
					try{
						forbiddenRegionImage = forbiddenRegionCartesian.getSelectedSequence();
						forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
					} catch (Exception eFRD) {
						forbiddenRegionArray = new int[0][0];
					}
				}
				if(InputTypeBox.getSelectedItem().toString()=="Distance to cell border"){
					CoordinateSystemLayout.show(CoordinateSystemPanel, "Distance");
					try{
						cellMaskDistance.getSelectedItem().toString();
						cellMaskImage = cellMaskDistance.getSelectedSequence();
						arraySize=cellMaskImage.getSizeX()*cellMaskImage.getSizeY();
						width=cellMaskImage.getSizeX();
						height=cellMaskImage.getSizeY();
						depth=maxZ+1;
						cellMaskArray = computeParametersAssociatedWithCellMask();
					} catch (Exception eFRD) {
						cellMaskArray = new int[0][0];
					}
					try{
						forbiddenRegionImage = forbiddenRegionDistance.getSelectedSequence();
						forbiddenRegionArray = computeParametersAssociatedWithForbiddenRegionMask();
					} catch (Exception eFRD) {
						forbiddenRegionArray = new int[0][0];
					}
				}
			}
		});
		String[] TrackFeaturePossibilities = {"Confinement ratio","Displacement distance","Total path length","Lifetime"};
		InputChoicePanel.add(new JLabel("Track feature"));
		TrackFeatureBox = new JComboBox<String>(TrackFeaturePossibilities);
		InputChoicePanel.add(TrackFeatureBox);
		northPanel.add(InputChoicePanel, BorderLayout.NORTH);
		northPanel.add(CoordinateSystemPanel, BorderLayout.CENTER);
		
		
		this.panel.add(northPanel, BorderLayout.NORTH);
		
		JPanel action = new JPanel(new FlowLayout(FlowLayout.CENTER));
		JButton export = new JButton("Export to Excel");
		export.addActionListener(new ActionListener()
		{
			@Override
			public void actionPerformed(ActionEvent e)
			{
				File file = getValidSaveFile();
				if (file == null)
	                return;
	            try
	            {
	            	// variable initialization
	            	float[] weight = new float[0],
	            			component1 = new float[0],
	            			component2 = new float[0],
	            			component3 = new float[0],
	            			distance1 = new float[0],
	            			distance2 = new float[0];
	            	int referenceCenter1,referenceCenter2,referenceCenter3,nbTrajectories=0;
	            	float referenceDirection1=(float)Math.PI/2;

	            	// outer circle definition from mask
	            	int[][] cellBorder = new int[0][0], innerBorder = new int[0][0];
	            	if(InputTypeBox.getSelectedItem().toString()=="Cylindrical"){
	            		try{
	            			cellMaskCylindrical.getSelectedItem().toString();
	            			cellBorder = computePseudo3DCellBorder(cellMaskArray,width,height,depth);
	            			depth=maxZ+1;
	            		} catch (Exception eFRD) {
	            			cellCenter1=(maxX-minX)/2;
	            			cellCenter2=(maxY-minY)/2;
	            			cellCenter3=(maxZ-minZ)/2;
	            			if(forbiddenRegionCylindrical.getSelectedSequence()==null){
	            				width=maxX+1;
	            				height=maxY+2;
	            				depth=maxZ+1;
	            				arraySize=width*height;
	            			}
	            			cellBorder = new int[depth][arraySize];
	            		}
	            		try{
	            			forbiddenRegionCylindrical.getSelectedItem().toString();
	            			innerBorder = computePseudo3DInnerBorder(forbiddenRegionArray,width,height,depth);
	            		} catch (Exception eFRD) {
	            			innerBorder = new int[depth][arraySize];
	            		}

	            		// reference center
	            		if(((int)(ref1Xcylindrical.getValue())>-1)&&((int)(ref1Xcylindrical.getValue())<width)&&((int)(ref1Ycylindrical.getValue())>-1)&&((int)(ref1Ycylindrical.getValue())<height)){
	            			referenceCenter1 = (int)(ref1Xcylindrical.getValue());
	            			referenceCenter2 = (int)(ref1Ycylindrical.getValue());
	            		}
	            		else{
	            			referenceCenter1 = cellCenter1;
	            			referenceCenter2 = cellCenter2;
	            		}

	            		// reference direction
	            		if(((int)(ref2Xcylindrical.getValue())>-1)&&((int)(ref2Xcylindrical.getValue())<width)&&((int)(ref2Ycylindrical.getValue())>-1)&&((int)(ref2Ycylindrical.getValue())<height)){
	            			referenceDirection1 = (float)(Math.atan2((float)(referenceCenter2-(int)(ref2Ycylindrical.getValue())),(float)((int)(ref2Xcylindrical.getValue())-referenceCenter1)) + Math.PI/2.);
	            		}
	            		else{
	            			if((referenceCenter1!=cellCenter1)&&(referenceCenter2!=cellCenter2)){
	            				referenceDirection1 = (float)(Math.atan2((float)(referenceCenter2-cellCenter2),(float)(cellCenter1-referenceCenter1)) + Math.PI/2.);
	            			}
	            		}
	            		
	            		// sum over time and extract intensity
	            		float[] tracksProjection = new float[arraySize];
	            		ArrayList<TrackGroup> TG1 = trackPool.getTrackGroupList();
	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG1.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG1.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			nbTrajectories += trackSegmentsA.size();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				tracksProjection[(int)(pt.getY())*width+(int)(pt.getX())] += 1;
	            			}
	            		}
	            		
	            		float maxRadius=0;
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
	            					if(tracksProjection[y*width+x]>0){
	            						if((cellMaskCylindrical.getSelectedSequence()!=null)||(forbiddenRegionCylindrical.getSelectedSequence()!=null)){
	            							if(forbiddenRegionCylindrical.getSelectedSequence()!=null){
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

	            		// extract tracks coordinates
	            		ArrayList<TrackGroup> TG2 = trackPool.getTrackGroupList();
	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG2.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG2.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				float currentTrackFeature = computeTrackFeature(tsA,TrackFeatureBox.getSelectedItem().toString());
	            				if(cellMaskCylindrical.getSelectedSequence()!=null){
	            					if(forbiddenRegionCylindrical.getSelectedSequence()!=null){
	            						if((cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0)&&(forbiddenRegionArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]==0)){
	            							// compute angle
	            							float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));
	            							// take into account the orientation reference
	            							float orientation=(theta-referenceDirection1);
	            							while(orientation<0){orientation+=(2*Math.PI);}
	            							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}
	            							// initialization
	            							float currentRadius=0;
	            							if(forbiddenRegionArray[(int)(pt.getZ())][referenceCenter2*width+referenceCenter1]>0){
	            								currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getX(),2.)+Math.pow(insideYCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getY(),2.)));
	            							}
	            							else{
	            								currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            							}
	            							// store coordinates 
	            							if(currentRadius>maxRadius){maxRadius = currentRadius;}
	            							extractCylindricalCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            					else{
	            						if(cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0){
	            							// compute angle
	            							float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));
	            							// take into account the orientation reference
	            							float orientation=(theta-referenceDirection1);
	            							while(orientation<0){orientation+=(2*Math.PI);}
	            							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}
	            							// initialization
	            							float currentRadius=0;
	            							currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            							// store coordinates 
	            							if(currentRadius>maxRadius){maxRadius = currentRadius;}
	            							extractCylindricalCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            				}
	            				else{
	            					// compute angle
	            					float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));

	            					// take into account the orientation reference
	            					float orientation=(theta-referenceDirection1);
	            					while(orientation<0){orientation+=(2*Math.PI);}
	            					while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            					// initialization
	            					float currentRadius=0;
	            					if(forbiddenRegionCylindrical.getSelectedSequence()!=null){
	            						if(forbiddenRegionArray[(int)(pt.getZ())][referenceCenter2*width+referenceCenter1]>0){
	            							currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getX(),2.)+Math.pow(insideYCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getY(),2.)));
	            						}
	            						else{
	            							currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            						}
	            					}
	            					else{
	            						currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            					}

	            					// store coordinates 
	            					if(currentRadius>maxRadius){maxRadius = currentRadius;}
	            					extractCylindricalCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            				}
	            			}
	            		}
	            		// radius
	            		// compute density for radius distribution for each experiment
	            		float bandwidthForGaussianDistributionForRadius = compute_rule_of_thumb_bandwidth(component1);

	            		// compute radius histograms
	            		// To modify distanceToZborder -> distanceToPlaneBorder
	            		float[] radiusHistogram = computeRadiusHistogram(component1,weight,distance1,(int)(nbBinsForCylindricalRadius.getValue()));

	            		// compute Gaussian kernel for each experiment
	            		float[] GaussianKernelForRadius = computeGaussianKernel((int)(nbBinsForCylindricalRadius.getValue()),bandwidthForGaussianDistributionForRadius,maxRadius);

	            		// compute radius density
	            		float[] radiusDensity = computeDensity((int)(nbBinsForCylindricalRadius.getValue()),radiusHistogram,GaussianKernelForRadius);

	            		// theta
	            		// compute density for theta distribution for each experiment
	            		float concentrationForVonMisesDistribution = compute_rule_of_thumb_concentration(component2);

	            		// compute theta histograms
	            		float[] angleHistogram = computeCircularHistogram(component2,weight,distance1,(int)(nbBinsForAngle.getValue()));

	            		// compute von Mises kernel for each experiment
	            		float[] vonMisesKernel = computeVonMisesKernel((int)(nbBinsForAngle.getValue()),concentrationForVonMisesDistribution);

	            		// compute circular density
	            		float[] angleDensity = computeCircularDensity((int)(nbBinsForAngle.getValue()),angleHistogram,vonMisesKernel);

	            		// depth
	            		// compute density for depth distribution for each experiment
	            		if(depth>1){
	            			// compute density for depth distribution for each experiment
	            			float bandwidthForGaussianDistributionForDepth = compute_rule_of_thumb_bandwidth(component3);

	            			// compute depth histograms
	            			float[] interpolatedDepthHistogram = new float[(int)(nbBinsForCylindricalDepth.getValue())],
	            					depthHistogram = computeDepthHistogram(component3,weight,distance2,(int)(nbBinsForCylindricalDepth.getValue()),depth-1,interpolatedDepthHistogram);

	            			// compute Gaussian kernel for each experiment
	            			float[] GaussianKernelForDepth = computeGaussianKernel((int)(nbBinsForCylindricalDepth.getValue()),bandwidthForGaussianDistributionForDepth,depth-1);

	            			// compute depth density
	            			float[] depthDensity = computeDensity((int)(nbBinsForCylindricalDepth.getValue()),depthHistogram,GaussianKernelForDepth);

	            			if (!FileUtil.getFileExtension(file.getPath(), false).equalsIgnoreCase("xls")) file = new File(file.getPath() + ".xls");
	            			exportCylindricalFiles(radiusHistogram,angleHistogram,depthHistogram,radiusDensity,angleDensity,depthDensity,file);

	            		}
	            		else{
	            			if (!FileUtil.getFileExtension(file.getPath(), false).equalsIgnoreCase("xls")) file = new File(file.getPath() + ".xls");
	            			exportCylindricalFiles(radiusHistogram,angleHistogram,radiusDensity,angleDensity,file);
	            		}
	            	}
	            	if(InputTypeBox.getSelectedItem().toString()=="Cartesian"){
	            		try{
	            			cellMaskCartesian.getSelectedItem().toString();
	            			cellBorder = computePseudo3DCellBorder(cellMaskArray,width,height,depth);
	            		} catch (Exception eFRD) {
	            			cellCenter1=(maxX-minX)/2;
	            			cellCenter2=(maxY-minY)/2;
	            			cellCenter3=(maxZ-minZ)/2;
	            			if(forbiddenRegionCartesian.getSelectedSequence()==null){
	            				width=maxX+1;
	            				height=maxY+2;
	            				depth=maxZ+1;
	            				arraySize=width*height;
	            			}
	            			cellBorder = new int[depth][arraySize];
	            		}

	            		try{
	            			forbiddenRegionCartesian.getSelectedItem().toString();
	            			innerBorder = computePseudo3DInnerBorder(forbiddenRegionArray,width,height,depth);
	            		} catch (Exception eFRD) {
	            			innerBorder = new int[depth][(maxX+1)*(maxY+1)];
	            		}
	            		// reference center
	            		if(((int)(ref1Xcartesian.getValue())>-1)&&((int)(ref1Xcartesian.getValue())<width)&&((int)(ref1Ycartesian.getValue())>-1)&&((int)(ref1Ycartesian.getValue())<height)){
	            			referenceCenter1 = (int)(ref1Xcartesian.getValue());
	            			referenceCenter2 = (int)(ref1Ycartesian.getValue());
	            		}
	            		else{
	            			referenceCenter1 = cellCenter1;
	            			referenceCenter2 = cellCenter2;
	            		}
	            		// reference direction
	            		if(((int)(ref2Xcartesian.getValue())>-1)&&((int)(ref2Xcartesian.getValue())<width)&&((int)(ref2Ycartesian.getValue())>-1)&&((int)(ref2Ycartesian.getValue())<height)){
	            			referenceDirection1 = (float)(Math.atan2((float)(referenceCenter2-(int)(ref2Ycartesian.getValue())),(float)((int)(ref2XcartesianSpinner.getValue())-referenceCenter1)) + Math.PI/2.);
	            		}
	            		else{
	            			if((referenceCenter1!=cellCenter1)&&(referenceCenter2!=cellCenter2)){
	            				referenceDirection1 = (float)(Math.atan2((float)(referenceCenter2-cellCenter2),(float)(cellCenter1-referenceCenter1)) + Math.PI/2.);
	            			}
	            		}
	            		// sum over time and extract intensity
	            		float[] tracksProjection = new float[arraySize];
	            		ArrayList<TrackGroup> TG = trackPool.getTrackGroupList();

	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			nbTrajectories += trackSegmentsA.size();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				tracksProjection[(int)(pt.getY())*width+(int)(pt.getX())] += 1;
	            			}
	            		}

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
	            					if(tracksProjection[y*width+x]>0.001){
	            						if((cellMaskCartesian.getSelectedSequence()!=null)||(forbiddenRegionCartesian.getSelectedSequence()!=null)){
	            							if(forbiddenRegionCartesian.getSelectedSequence()!=null){
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
	            		// extract tracks coordinates
	            		ArrayList<TrackGroup> TG2 = trackPool.getTrackGroupList();
	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG2.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG2.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				float currentTrackFeature = computeTrackFeature(tsA,TrackFeatureBox.getSelectedItem().toString());
	            				if(cellMaskCartesian.getSelectedSequence()!=null){
	            					if(forbiddenRegionCartesian.getSelectedSequence()!=null){
	            						if((cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0)&&(forbiddenRegionArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]==0)){
	            							// compute angle
	            							float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));

	            							// take into account the orientation reference
	            							float orientation=(theta-referenceDirection1);
	            							while(orientation<0){orientation+=(2*Math.PI);}
	            							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            							// initialization
	            							float currentRadius=0;
	            							if(forbiddenRegionArray[(int)(pt.getZ())][referenceCenter2*width+referenceCenter1]>0){
	            								currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getX(),2.)+Math.pow(insideYCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getY(),2.)));
	            							}
	            							else{
	            								currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            							}
	            							// store coordinates 
	            							extractCartesianCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            					else{
	            						if(cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0){
	            							// compute angle
	            							float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));

	            							// take into account the orientation reference
	            							float orientation=(theta-referenceDirection1);
	            							while(orientation<0){orientation+=(2*Math.PI);}
	            							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            							// initialization
	            							float currentRadius=0;
	            							currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            							// store coordinates 
	            							extractCartesianCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            				}
	            				else{
	            					// compute angle
	            					float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));

	            					// take into account the orientation reference
	            					float orientation=(theta-referenceDirection1);
	            					while(orientation<0){orientation+=(2*Math.PI);}
	            					while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            					// initialization
	            					float currentRadius=0;
	            					if(forbiddenRegionCartesian.getSelectedSequence()!=null){
	            						if(forbiddenRegionArray[(int)(pt.getZ())][referenceCenter2*width+referenceCenter1]>0){
	            							currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getX(),2.)+Math.pow(insideYCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getY(),2.)));
	            						}
	            						else{
	            							currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            						}
	            					}
	            					else{
	            						currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            					}

	            					// store coordinates 
	            					extractCartesianCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            				}
	            			}
	            		}
	            		// x
	            		// compute density for x distribution for each experiment
	            		float bandwidthForGaussianDistributionForX = compute_rule_of_thumb_bandwidth(component1);

	            		// compute x histograms
	            		float[] interpolatedXHistogram = new float[(int)(nbBinsForX.getValue())],
	            				xHistogram = computeSymmetricalHistogram(component1,weight,distance1,(int)(nbBinsForX.getValue()),interpolatedXHistogram);

	            		// compute Gaussian kernel
	            		float[] GaussianKernelForX = computeSymmetricalGaussianKernel((int)(nbBinsForX.getValue()),bandwidthForGaussianDistributionForX,component1);

	            		// compute x density
	            		float[] xDensity = computeDensity((int)(nbBinsForX.getValue()),xHistogram,GaussianKernelForX);

	            		// y
	            		// compute density for y distribution for each experiment
	            		float bandwidthForGaussianDistributionForY = compute_rule_of_thumb_bandwidth(component2);

	            		// compute x histograms
	            		float[] interpolatedYHistogram = new float[(int)(nbBinsForY.getValue())],
	            				yHistogram = computeSymmetricalHistogram(component2,weight,distance1,(int)(nbBinsForY.getValue()),interpolatedYHistogram);

	            		// compute Gaussian kernel
	            		float[] GaussianKernelForY = computeSymmetricalGaussianKernel((int)(nbBinsForY.getValue()),bandwidthForGaussianDistributionForY,component2);

	            		// compute x density
	            		float[] yDensity = computeDensity((int)(nbBinsForY.getValue()),yHistogram,GaussianKernelForY);

	            		// depth
	            		// compute density for depth distribution for each experiment
	            		if(depth>1){
	            			// compute density for depth distribution for each experiment
	            			float bandwidthForGaussianDistributionForDepth = compute_rule_of_thumb_bandwidth(component3);

	            			// compute depth histograms
	            			float[] interpolatedDepthHistogram = new float[(int)(nbBinsForCartesianDepth.getValue())],
	            					depthHistogram = computeDepthHistogram(component3,weight,distance2,(int)(nbBinsForCartesianDepth.getValue()),depth-1,interpolatedDepthHistogram);

	            			// compute Gaussian kernel for each experiment
	            			float[] GaussianKernelForDepth = computeGaussianKernel((int)(nbBinsForCartesianDepth.getValue()),bandwidthForGaussianDistributionForDepth,depth-1);

	            			// compute depth density
	            			float[] depthDensity = computeDensity((int)(nbBinsForCartesianDepth.getValue()),depthHistogram,GaussianKernelForDepth);

	            			if (!FileUtil.getFileExtension(file.getPath(), false).equalsIgnoreCase("xls")) file = new File(file.getPath() + ".xls");
	            			exportCartesianFiles(xHistogram,yHistogram,depthHistogram,xDensity,yDensity,depthDensity,file);
	            		}
	            		else{
	            			if (!FileUtil.getFileExtension(file.getPath(), false).equalsIgnoreCase("xls")) file = new File(file.getPath() + ".xls");
	            			exportCartesianFiles(xHistogram,yHistogram,xDensity,yDensity,file);
	            		}
	            	}
	            	if(InputTypeBox.getSelectedItem().toString()=="Spherical"){
	            		try{
	            			cellMaskSpherical.getSelectedItem().toString();
	            			cellBorder = computePseudo3DCellBorder(cellMaskArray,width,height,depth);
	            		} catch (Exception eFRD) {
	            			cellBorder = new int[1][(maxX+1)*(maxY+1)];
	            			if(forbiddenRegionSpherical.getSelectedSequence()==null){
	            				width=maxX+1;
	            				height=maxY+2;
	            				depth=maxZ+1;
	            				arraySize=width*height;
	            			}
	            			cellBorder = new int[depth][arraySize];
	            		}

	            		try{
	            			forbiddenRegionSpherical.getSelectedItem().toString();
	            			innerBorder = computePseudo3DInnerBorder(forbiddenRegionArray,width,height,depth);
	            		} catch (Exception eFRD) {
	            			innerBorder = new int[depth][(maxX+1)*(maxY+1)];
	            		}
	            		// reference center
	            		if(((int)(ref1Xspherical.getValue())>-1)&&((int)(ref1Xspherical.getValue())<width)&&((int)(ref1Yspherical.getValue())>-1)&&((int)(ref1Yspherical.getValue())<height)){
	            			referenceCenter1 = (int)(ref1Xspherical.getValue());
	            			referenceCenter2 = (int)(ref1Yspherical.getValue());
	            			if(((int)(ref1Zspherical.getValue())>-1)&&((int)(ref1Zspherical.getValue())<depth)){
	            				referenceCenter3 = (int)(ref1Zspherical.getValue());
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
	            		if(((int)(ref2Xspherical.getValue())>-1)&&((int)(ref2Xspherical.getValue())<width)&&((int)(ref2Yspherical.getValue())>-1)&&((int)(ref2Yspherical.getValue())<height)){
	            			referenceDirection1 = (float)(Math.atan2((referenceCenter2-(int)(ref2Yspherical.getValue())),((int)(ref2Xspherical.getValue())-referenceCenter1)) + Math.PI/2.);
	            		}
	            		else{
	            			if((referenceCenter1!=cellCenter1)&&(referenceCenter2!=cellCenter2)){
	            				referenceDirection1 = (float)(Math.atan2((referenceCenter2-cellCenter2),(cellCenter1-referenceCenter1)) + Math.PI/2.);
	            			}
	            		}
	            		// sum over time and extract intensity
	            		float[] tracksProjection = new float[arraySize];
	            		ArrayList<TrackGroup> TG = trackPool.getTrackGroupList();

	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			nbTrajectories += trackSegmentsA.size();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				tracksProjection[(int)(pt.getY())*width+(int)(pt.getX())] += 1;
	            			}
	            		}
	            		float maxRadius=0;
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
	            					if(tracksProjection[y*width+x]>0.001){
	            						if((cellMaskSpherical.getSelectedSequence()!=null)||(forbiddenRegionSpherical.getSelectedSequence()!=null)){
	            							if(forbiddenRegionSpherical.getSelectedSequence()!=null){
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
	            		
	            		// extract tracks coordinates
	            		ArrayList<TrackGroup> TG2 = trackPool.getTrackGroupList();
	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG2.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG2.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				float currentTrackFeature = computeTrackFeature(tsA,TrackFeatureBox.getSelectedItem().toString());
	            				if(cellMaskSpherical.getSelectedSequence()!=null){
	            					if(forbiddenRegionSpherical.getSelectedSequence()!=null){
	            						if((cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0)&&(forbiddenRegionArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]==0)){
	            							// compute angle
	            							float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));
	            							// take into account the orientation reference
	            							float orientation=(theta-referenceDirection1);
	            							while(orientation<0){orientation+=(2*Math.PI);}
	            							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            							// initialization
	            							float currentRadius=0;
	            							if(forbiddenRegionArray[(int)(pt.getZ())][referenceCenter2*width+referenceCenter1]>0){
	            								currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getX(),2.)+Math.pow(insideYCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getY(),2.)));
	            							}
	            							else{
	            								currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            							}
	            							// store coordinates 
	            							if(currentRadius>maxRadius){maxRadius = currentRadius;}
	            							extractSphericalCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            					else{
	            						if(cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0){
	            							// compute angle
	            							float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));

	            							// take into account the orientation reference
	            							float orientation=(theta-referenceDirection1);
	            							while(orientation<0){orientation+=(2*Math.PI);}
	            							while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            							// initialization
	            							float currentRadius=0;
	            							currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            							// store coordinates 
	            							if(currentRadius>maxRadius){maxRadius = currentRadius;}
	            							extractSphericalCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            				}
	            				else{
	            					// compute angle
	            					float theta=(float)(Math.atan2((float)(referenceCenter2)-pt.getY(),pt.getX()-(float)(referenceCenter1)));

	            					// take into account the orientation reference
	            					float orientation=(theta-referenceDirection1);
	            					while(orientation<0){orientation+=(2*Math.PI);}
	            					while(orientation>(2*Math.PI)){orientation-=(2*Math.PI);}

	            					// initialization
	            					float currentRadius=0;
	            					if(forbiddenRegionSpherical.getSelectedSequence()!=null){
	            						if(forbiddenRegionArray[(int)(pt.getZ())][referenceCenter2*width+referenceCenter1]>0){
	            							currentRadius =	(float)(Math.sqrt(Math.pow(insideXCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getX(),2.)+Math.pow(insideYCoord[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]-pt.getY(),2.)));
	            						}
	            						else{
	            							currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            						}
	            					}
	            					else{
	            						currentRadius =	(float)(Math.sqrt(Math.pow((float)(referenceCenter1)-pt.getX(),2.)+Math.pow((float)(referenceCenter2)-pt.getY(),2.)));
	            					}

	            					// store coordinates 
	            					if(currentRadius>maxRadius){maxRadius = currentRadius;}
	            					extractSphericalCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), currentRadius, orientation, cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, component2, component3, weight, distance1, width, height, depth, trackSegmentIndexA);
	            				}
	            			}
	            		}
	            		// radius
	            		// compute density for radius distribution for each experiment
	            		float bandwidthForGaussianDistributionForRadius = compute_rule_of_thumb_bandwidth(component1);

	            		// compute radius histograms
	            		// To modify distanceToZborder -> distanceToPlaneBorder
	            		float[] radiusHistogram = computeRadiusHistogram(component1,weight,distance1,(int)(nbBinsForSphericalRadius.getValue()));

	            		// compute Gaussian kernel for each experiment
	            		float[] GaussianKernelForRadius = computeGaussianKernel((int)(nbBinsForSphericalRadius.getValue()),bandwidthForGaussianDistributionForRadius,maxRadius);

	            		// compute radius density
	            		float[] radiusDensity = computeDensity((int)(nbBinsForSphericalRadius.getValue()),radiusHistogram,GaussianKernelForRadius);


	            		// azimuthal angle
	            		// compute density for azimuthal angle distribution for each experiment
	            		float concentrationForVonMisesDistributionOfAzimuthalAngle = compute_rule_of_thumb_concentration(component2);

	            		// compute azimuthal angle histograms
	            		float[] azimuthalAngleHistogram = computeCircularHistogram(component2,weight,distance1,(int)(nbBinsForFirstAngle.getValue()));

	            		// compute von Mises kernel for each experiment
	            		float[] vonMisesKernelOfAzimuthalAngle = computeVonMisesKernel((int)(nbBinsForFirstAngle.getValue()),concentrationForVonMisesDistributionOfAzimuthalAngle);

	            		// compute azimuthal angle density
	            		float[] azimuthalAngleDensity = computeCircularDensity((int)(nbBinsForFirstAngle.getValue()),azimuthalAngleHistogram,vonMisesKernelOfAzimuthalAngle);


	            		// polar angle
	            		// compute density for polar angle distribution for each experiment
	            		float concentrationForVonMisesDistributionOfPolarAngle = compute_rule_of_thumb_concentration(component3);

	            		// compute polar angle histograms
	            		float[] polarAngleHistogram = computeCircularHistogram(component3,weight,distance1,(int)(nbBinsForSecondAngle.getValue()));

	            		// compute von Mises kernel for each experiment
	            		float[] vonMisesKernelOfPolarAngle = computeVonMisesKernel((int)(nbBinsForSecondAngle.getValue()),concentrationForVonMisesDistributionOfPolarAngle);

	            		// compute polar angle density
	            		float[] polarAngleDensity = computeCircularDensity((int)(nbBinsForSecondAngle.getValue()),polarAngleHistogram,vonMisesKernelOfPolarAngle);

	            		if (!FileUtil.getFileExtension(file.getPath(), false).equalsIgnoreCase("xls")) file = new File(file.getPath() + ".xls");
	            		exportSphericalFiles(radiusHistogram,azimuthalAngleHistogram,polarAngleHistogram,radiusDensity,azimuthalAngleDensity,polarAngleDensity,file);
	            	}
	            	if(InputTypeBox.getSelectedItem().toString()=="Distance to cell border"){
	            		try{
	            			cellMaskDistance.getSelectedItem().toString();
	            			cellBorder = computePseudo3DCellBorder(cellMaskArray,width,height,depth);
	            			depth=maxZ+1;
	            		} catch (Exception eFRD) {
	            			cellCenter1=(maxX-minX)/2;
	            			cellCenter2=(maxY-minY)/2;
	            			cellCenter3=(maxZ-minZ)/2;
	            			if(forbiddenRegionDistance.getSelectedSequence()==null){
	            				width=maxX+1;
	            				height=maxY+2;
	            				depth=maxZ+1;
	            				arraySize=width*height;
	            			}
	            			cellBorder = new int[depth][arraySize];
	            		}

	            		try{
	            			forbiddenRegionDistance.getSelectedItem().toString();
	            			innerBorder = computePseudo3DInnerBorder(forbiddenRegionArray,width,height,depth);
	            		} catch (Exception eFRD) {
	            			innerBorder = new int[depth][(maxX+1)*(maxY+1)];
	            		}
	            		// sum over time and extract intensity
	            		float[] tracksProjection = new float[arraySize];
	            		ArrayList<TrackGroup> TG = trackPool.getTrackGroupList();

	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			nbTrajectories += trackSegmentsA.size();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				tracksProjection[(int)(pt.getY())*width+(int)(pt.getX())] += 1;
	            			}
	            		}

	            		float maxDistance=0;
	            		component1 = new float[nbTrajectories];
	            		weight = new float[nbTrajectories];
	            		distance1 = new float[nbTrajectories];

	            		// compute normalizing distance
	            		float[][] cellSegment = new float[depth][arraySize];
	            		for(int z=0;z<depth;z++){
	            			for(int y=0;y<height;y++){
	            				for(int x=0;x<width;x++){
	            					if(tracksProjection[y*width+x]>0.001){
	            						if((cellMaskDistance.getSelectedSequence()!=null)||(forbiddenRegionDistance.getSelectedSequence()!=null)){
	            							if(forbiddenRegionDistance.getSelectedSequence()!=null){
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
	            		// extract tracks coordinates
	            		ArrayList<TrackGroup> TG2 = trackPool.getTrackGroupList();
	            		for (int trackGroupIndexA = 0; trackGroupIndexA < TG2.size(); trackGroupIndexA++)
	            		{
	            			TrackGroup tgA = TG2.get(trackGroupIndexA);
	            			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();
	            			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
	            			{
	            				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
	            				Detection pt = tsA.getDetectionList().get(tsA.getDetectionList().size()/2);
	            				float currentTrackFeature = computeTrackFeature(tsA,TrackFeatureBox.getSelectedItem().toString());
	            				if(cellMaskDistance.getSelectedSequence()!=null){
	            					if(forbiddenRegionDistance.getSelectedSequence()!=null){
	            						if((cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0)&&(forbiddenRegionArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]==0)){
	            							extractDistanceCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, weight, distance1, width, height, depth, trackSegmentIndexA);	
	            						}
	            					}
	            					else{
	            						if(cellMaskArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>0){
	            							extractDistanceCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            				}
	            				else{
	            					if(forbiddenRegionDistance.getSelectedSequence()!=null){
	            						if(forbiddenRegionArray[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]==0){
	            							extractDistanceCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, weight, distance1, width, height, depth, trackSegmentIndexA);
	            						}
	            					}
	            					else{
	            						extractDistanceCoordinates((int)(pt.getX()), (int)(pt.getY()), (int)(pt.getZ()), cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())], currentTrackFeature, component1, weight, distance1, width, height, depth, trackSegmentIndexA);
	            					}
	            				}
	            				if(cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())]>maxDistance){maxDistance = cellSegment[(int)(pt.getZ())][(int)(pt.getY())*width+(int)(pt.getX())];}
	            			}
	            		}
	            		// distance to cell border
	            		// compute density for radius distribution for each experiment
	            		float bandwidthForGaussianDistributionForDistance = compute_rule_of_thumb_bandwidth(component1);

	            		// compute distance histograms
	            		float[] distanceHistogram = computeRadiusHistogram(component1,weight,distance1,(int)(nbBinsForDistance.getValue()));

	            		// compute Gaussian kernel for each experiment
	            		float[] GaussianKernelForDistance = computeGaussianKernel((int)(nbBinsForDistance.getValue()),bandwidthForGaussianDistributionForDistance,maxDistance);

	            		// compute distance density
	            		float[] distanceDensity = computeDensity((int)(nbBinsForDistance.getValue()),distanceHistogram,GaussianKernelForDistance);

	            		if (!FileUtil.getFileExtension(file.getPath(), false).equalsIgnoreCase("xls")) file = new File(file.getPath() + ".xls");
	            		exportDistanceFiles(distanceHistogram,distanceDensity,file);
	            	}
	            }catch (Exception eFile)
	            {
	            	MessageDialog.showDialog("Writing the save file failed.");
	            }
			}
		});
		action.add(export);
		this.panel.add(action, BorderLayout.SOUTH);
		this.panel.setVisible(true);
		
		ArrayList<TrackGroup> trackGroups = trackPool.getTrackGroupList();

		for (int trackGroupIndexA = 0; trackGroupIndexA < trackGroups.size(); trackGroupIndexA++)
		{
			TrackGroup tgA = trackGroups.get(trackGroupIndexA);
			ArrayList<TrackSegment> trackSegmentsA = tgA.getTrackSegmentList();

			for (int trackSegmentIndexA = 0; trackSegmentIndexA < trackSegmentsA.size(); trackSegmentIndexA++)
			{
				TrackSegment tsA = trackSegmentsA.get(trackSegmentIndexA);
				for (Detection dA : tsA.getDetectionList())
				{
					if((int)(dA.getX())<minX){minX=(int)(dA.getX());}
					if((int)(dA.getX())>maxX){maxX=(int)(dA.getX());}
					if((int)(dA.getY())<minY){minY=(int)(dA.getY());}
					if((int)(dA.getY())>maxY){maxY=(int)(dA.getY());}
					if((int)(dA.getZ())<minZ){minZ=(int)(dA.getZ());}
					if((int)(dA.getZ())>maxZ){maxZ=(int)(dA.getZ());}
				}
			}
		}
	}

	@Override
	public void displaySequenceChanged()
	{
	}

}
