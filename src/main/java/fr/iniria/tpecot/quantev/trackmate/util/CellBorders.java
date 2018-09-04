package fr.iniria.tpecot.quantev.trackmate.util;

import static fr.iniria.tpecot.quantev.trackmate.util.Morphology.dilate;
import static fr.iniria.tpecot.quantev.trackmate.util.Morphology.makeStrElt;

public class CellBorders
{

	/**
	 * Compute cell border.
	 *
	 * @param mask
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @return
	 */
	public final static int[][] computePseudo3DCellBorder( final int[][] mask, final int dimX, final int dimY, final int dimZ )
	{
		final int strEltXY = 2;
		int strEltZ = 1;
		if ( dimZ > 1 )
		{
			strEltZ = 2;
		}
		final int[][] strElt = makeStrElt( strEltXY, strEltZ );
		final int[][] output = dilate( mask, strElt, strEltXY, strEltZ, dimX, dimY, dimZ );
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int i = 0; i < ( dimX * dimY ); i++ )
			{
				output[ z ][ i ] = output[ z ][ i ] - mask[ z ][ i ];
			}
		}

		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				if ( mask[ z ][ y * dimX ] == 1 )
				{
					output[ z ][ y * dimX ] = 1;
				}
				if ( mask[ z ][ y * dimX + dimX - 1 ] == 1 )
				{
					output[ z ][ y * dimX + dimX - 1 ] = 1;
				}
			}
			for ( int x = 0; x < dimX; x++ )
			{
				if ( mask[ z ][ x ] == 1 )
				{
					output[ z ][ x ] = 1;
				}
				if ( mask[ z ][ ( dimY - 1 ) * dimX + x ] == 1 )
				{
					output[ z ][ ( dimY - 1 ) * dimX + x ] = 1;
				}
			}
		}
		return output;
	}

	public static final int[] extractXcellBorderCoordinates( final int[][] cellBorder, final int dimX, final int dimY, final int dimZ )
	{
		int nbTrajectories = 0;
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				for ( int x = 0; x < dimX; x++ )
				{
					if ( cellBorder[ z ][ y * dimX + x ] > 0 )
					{
						nbTrajectories++;
					}
				}
			}
		}
		final int[] output = new int[ nbTrajectories ];
		int i = 0;
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				for ( int x = 0; x < dimX; x++ )
				{
					if ( cellBorder[ z ][ y * dimX + x ] > 0 )
					{
						output[ i ] = x;
						i++;
					}
				}
			}
		}
		return output;
	}

	public static final int[] extractYcellBorderCoordinates( final int[][] cellBorder, final int dimX, final int dimY, final int dimZ )
	{
		int nbTrajectories = 0;
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				for ( int x = 0; x < dimX; x++ )
				{
					if ( cellBorder[ z ][ y * dimX + x ] > 0 )
					{
						nbTrajectories++;
					}
				}
			}
		}
		final int[] output = new int[ nbTrajectories ];
		int i = 0;
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				for ( int x = 0; x < dimX; x++ )
				{
					if ( cellBorder[ z ][ y * dimX + x ] > 0 )
					{
						output[ i ] = y;
						i++;
					}
				}
			}
		}
		return output;
	}

	/**
	 * Compute cell border.
	 *
	 * @param mask
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @return
	 */
	public static final int[][] compute3DCellBorder( final int[][] mask, final int dimX, final int dimY, final int dimZ )
	{
		final int strEltXY = 2;
		int strEltZ = 1;
		if ( dimZ > 1 )
		{
			strEltZ = 2;
		}
		final int[][] strElt = makeStrElt( strEltXY, strEltZ );
		final int[][] output = dilate( mask, strElt, strEltXY, strEltZ, dimX, dimY, dimZ );
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int i = 0; i < ( dimX * dimY ); i++ )
			{
				output[ z ][ i ] = output[ z ][ i ] - mask[ z ][ i ];
			}
		}

		for ( int i = 0; i < ( dimX * dimY ); i++ )
		{
			if ( mask[ 0 ][ i ] == 1 )
			{
				output[ 0 ][ i ] = 1;
			}
			if ( mask[ dimZ - 1 ][ i ] == 1 )
			{
				output[ dimZ - 1 ][ i ] = 1;
			}
		}
		for ( int z = 1; z < ( dimZ - 1 ); z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				if ( mask[ z ][ y * dimX ] == 1 )
				{
					output[ z ][ y * dimX ] = 1;
				}
				if ( mask[ z ][ y * dimX + dimX - 1 ] == 1 )
				{
					output[ z ][ y * dimX + dimX - 1 ] = 1;
				}
			}
			for ( int x = 0; x < dimX; x++ )
			{
				if ( mask[ z ][ x ] == 1 )
				{
					output[ z ][ x ] = 1;
				}
				if ( mask[ z ][ ( dimY - 1 ) * dimX + x ] == 1 )
				{
					output[ z ][ ( dimY - 1 ) * dimX + x ] = 1;
				}
			}
		}

		return output;
	}

	public static final int[][] computePseudo3DInnerBorder( final int[][] mask, final int dimX, final int dimY, final int dimZ )
	{
		final int strEltXY = 2;
		int strEltZ = 1;
		if ( dimZ > 1 )
		{
			strEltZ = 2;
		}
		final int[][] strElt = makeStrElt( strEltXY, strEltZ );
		final int[][] output = dilate( mask, strElt, strEltXY, strEltZ, dimX, dimY, dimZ );
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int i = 0; i < ( dimX * dimY ); i++ )
			{
				output[ z ][ i ] = output[ z ][ i ] - mask[ z ][ i ];
			}
		}
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				if ( mask[ z ][ y * dimX ] == 1 )
				{
					output[ z ][ y * dimX ] = 1;
				}
				if ( mask[ z ][ y * dimX + dimX - 1 ] == 1 )
				{
					output[ z ][ y * dimX + dimX - 1 ] = 1;
				}
			}
			for ( int x = 0; x < dimX; x++ )
			{
				if ( mask[ z ][ x ] == 1 )
				{
					output[ z ][ x ] = 1;
				}
				if ( mask[ z ][ ( dimY - 1 ) * dimX + x ] == 1 )
				{
					output[ z ][ ( dimY - 1 ) * dimX + x ] = 1;
				}
			}
		}

		return output;
	}

	public static final int[][] compute3DInnerBorder( final int[][] mask, final int dimX, final int dimY, final int dimZ )
	{
		final int strEltXY = 2;
		int strEltZ = 1;
		if ( dimZ > 1 )
		{
			strEltZ = 2;
		}
		final int[][] strElt = makeStrElt( strEltXY, strEltZ );
		final int[][] output = dilate( mask, strElt, strEltXY, strEltZ, dimX, dimY, dimZ );
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int i = 0; i < ( dimX * dimY ); i++ )
			{
				output[ z ][ i ] = output[ z ][ i ] - mask[ z ][ i ];
			}
		}
		for ( int i = 0; i < ( dimX * dimY ); i++ )
		{
			if ( mask[ 0 ][ i ] == 1 )
			{
				output[ 0 ][ i ] = 1;
			}
			if ( mask[ dimZ - 1 ][ i ] == 1 )
			{
				output[ dimZ - 1 ][ i ] = 1;
			}
		}
		for ( int z = 1; z < ( dimZ - 1 ); z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				if ( mask[ z ][ y * dimX ] == 1 )
				{
					output[ z ][ y * dimX ] = 1;
				}
				if ( mask[ z ][ y * dimX + dimX - 1 ] == 1 )
				{
					output[ z ][ y * dimX + dimX - 1 ] = 1;
				}
			}
			for ( int x = 0; x < dimX; x++ )
			{
				if ( mask[ z ][ x ] == 1 )
				{
					output[ z ][ x ] = 1;
				}
				if ( mask[ z ][ ( dimY - 1 ) * dimX + x ] == 1 )
				{
					output[ z ][ ( dimY - 1 ) * dimX + x ] = 1;
				}
			}
		}

		return output;
	}

	/**
	 * Compute distance to cell border.
	 *
	 * @param referenceCenter1
	 * @param referenceCenter2
	 * @param x
	 * @param y
	 * @param z
	 * @param cellBorder
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @return
	 */
	public static final float computeDistanceToCellBorder( final int referenceCenter1, final int referenceCenter2, final int x, final int y, final int z, final int[][] cellBorder, final int dimX, final int dimY, final int dimZ )
	{
		boolean verticalLine = false;
		float a = 0, b = 0;
		if ( x != referenceCenter1 )
		{
			a = ( ( float ) ( y - referenceCenter2 ) / ( float ) ( referenceCenter1 - x ) );
			b = ( -referenceCenter2 - a * referenceCenter1 );
		}
		else
		{
			verticalLine = true;
		}

		float minOutsideDistance = 10000;
		int outsideIdX = 0, outsideIdY = 0;

		for ( int j = 0; j < dimY; j++ )
		{
			for ( int i = 0; i < dimX; i++ )
			{
				if ( cellBorder[ z ][ j * dimX + i ] > 0 )
				{
					if ( ( Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) ) + Math.sqrt( Math.pow( x - referenceCenter1, 2. ) + Math.pow( y - referenceCenter2, 2. ) ) ) < ( Math.sqrt( Math.pow( i - referenceCenter1, 2. ) + Math.pow( j - referenceCenter2, 2. ) ) ) + 1. )
					{
						float distance = 0;
						if ( verticalLine )
						{
							if ( i == x )
							{
								distance = ( float ) ( Math.sqrt( Math.pow( j - referenceCenter2, 2. ) ) );
								if ( distance < minOutsideDistance )
								{
									minOutsideDistance = distance;
									outsideIdX = i;
									outsideIdY = j;
								}
							}
						}
						else
						{
							distance = ( float ) ( Math.abs( a * i + j + b ) / Math.sqrt( 1 + a * a ) );
							if ( distance < minOutsideDistance )
							{
								minOutsideDistance = distance;
								outsideIdX = i;
								outsideIdY = j;
							}
						}
					}
				}
			}
		}
		return ( float ) ( Math.sqrt( Math.pow( referenceCenter1 - outsideIdX, 2. ) + Math.pow( referenceCenter2 - outsideIdY, 2. ) ) );
	}

	public static final float computeDistanceToCellBorder( final int referenceCenter1, final int referenceCenter2, final int x, final int y, final int z, final int[][] cellBorder, final int[][] innerBorder, final int dimX, final int dimY, final int dimZ, final int[][] insideCoordX, final int[][] insideCoordY )
	{
		boolean verticalLine = false;
		float a = 0, b = 0;
		if ( x != referenceCenter1 )
		{
			a = ( ( float ) ( y - referenceCenter2 ) / ( float ) ( referenceCenter1 - x ) );
			b = ( -referenceCenter2 - a * referenceCenter1 );
		}
		else
		{
			verticalLine = true;
		}

		float minOutsideDistance = 10000, minInsideDistance = 10000;
		int outsideIdX = 0, outsideIdY = 0, insideIdX = 0, insideIdY = 0;

		for ( int j = 0; j < dimY; j++ )
		{
			for ( int i = 0; i < dimX; i++ )
			{
				if ( innerBorder[ z ][ j * dimX + i ] > 0 )
				{
					if ( ( Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) ) + Math.sqrt( Math.pow( i - referenceCenter1, 2. ) + Math.pow( j - referenceCenter2, 2. ) ) ) < ( Math.sqrt( Math.pow( x - referenceCenter1, 2. ) + Math.pow( y - referenceCenter2, 2. ) ) ) + 1. )
					{
						float distance = 0;
						if ( verticalLine )
						{
							if ( i == x )
							{
								distance = ( float ) ( Math.sqrt( Math.pow( j - referenceCenter2, 2. ) ) );
							}
						}
						else
						{
							distance = ( float ) ( Math.abs( a * i + j + b ) / Math.sqrt( 1 + a * a ) );
						}
						if ( distance < minInsideDistance )
						{
							minInsideDistance = distance;
							insideIdX = i;
							insideIdY = j;
						}
					}
				}
				if ( cellBorder[ z ][ j * dimX + i ] > 0 )
				{
					if ( ( Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) ) + Math.sqrt( Math.pow( x - referenceCenter1, 2. ) + Math.pow( y - referenceCenter2, 2. ) ) ) < ( Math.sqrt( Math.pow( i - referenceCenter1, 2. ) + Math.pow( j - referenceCenter2, 2. ) ) ) + 1. )
					{
						float distance = 0;
						if ( verticalLine )
						{
							if ( i == x )
							{
								distance = ( float ) ( Math.sqrt( Math.pow( j - referenceCenter2, 2. ) ) );
							}
						}
						else
						{
							distance = ( float ) ( Math.abs( a * i + j + b ) / Math.sqrt( 1 + a * a ) );
						}
						if ( distance < minOutsideDistance )
						{
							minOutsideDistance = distance;
							outsideIdX = i;
							outsideIdY = j;
						}

					}
				}
			}
		}
		insideCoordX[ z ][ y * dimX + x ] = insideIdX;
		insideCoordY[ z ][ y * dimX + x ] = insideIdY;
		return ( float ) ( Math.sqrt( Math.pow( insideIdX - outsideIdX, 2. ) + Math.pow( insideIdY - outsideIdY, 2. ) ) );
	}

	public static final float[] computeDistanceToCellBorder( final int[] cellBorderX, final int[] cellBorderY, final int nbTrajectories, final int x, final int y )
	{
		final float[] output = new float[ nbTrajectories * 3 ];
		final float angleUnit = ( float ) ( 2 * Math.PI / nbTrajectories );
		for ( int u = 0; u < cellBorderX.length; u++ )
		{
			float theta = ( float ) Math.atan2( y - cellBorderY[ u ], cellBorderX[ u ] - x );
			if ( theta < 0 )
			{
				theta += ( 2 * Math.PI );
			}
			if ( theta > ( 2 * Math.PI ) )
			{
				theta -= ( 2 * Math.PI );
			}
			output[ ( int ) ( theta / angleUnit ) ] = theta;
			output[ ( int ) ( theta / angleUnit ) + nbTrajectories ] = cellBorderX[ u ];
			output[ ( int ) ( theta / angleUnit ) + 2 * nbTrajectories ] = cellBorderY[ u ];
		}
		for ( int i = 0; i < nbTrajectories; i++ )
		{
			if ( ( output[ i ] < 0.001 ) && ( output[ i + nbTrajectories ] < 0.001 ) && ( output[ i + 2 * nbTrajectories ] < 0.001 ) )
			{
				int previousIndex = -1, nextIndex = -1;
				boolean previousOver = false, nextOver = false;
				for ( int u = i; u >= 0; u-- )
				{
					if ( ( output[ u + nbTrajectories ] > 0.001 ) && ( output[ u + 2 * nbTrajectories ] > 0.001 ) && ( !previousOver ) )
					{
						previousIndex = u;
						previousOver = true;
					}
				}
				for ( int u = i; u < nbTrajectories; u++ )
				{
					if ( ( output[ u + nbTrajectories ] > 0.001 ) && ( output[ u + 2 * nbTrajectories ] > 0.001 ) && ( !nextOver ) )
					{
						nextIndex = u;
						nextOver = true;
					}
				}
				output[ i ] = angleUnit * i;
				if ( ( previousIndex > -1 ) && ( nextIndex > -1 ) )
				{
					output[ i + nbTrajectories ] = ( output[ previousIndex + nbTrajectories ] + output[ nextIndex + nbTrajectories ] ) / 2;
					output[ i + 2 * nbTrajectories ] = ( output[ previousIndex + 2 * nbTrajectories ] + output[ nextIndex + 2 * nbTrajectories ] ) / 2;
				}
				else
				{
					if ( previousIndex > -1 )
					{
						output[ i + nbTrajectories ] = output[ previousIndex + nbTrajectories ];
						output[ i + 2 * nbTrajectories ] = output[ previousIndex + 2 * nbTrajectories ];
					}
					else
					{
						if ( nextIndex > -1 )
						{
							output[ i + nbTrajectories ] = output[ nextIndex + nbTrajectories ];
							output[ i + 2 * nbTrajectories ] = output[ nextIndex + 2 * nbTrajectories ];
						}
					}
				}
			}
		}
		return output;
	}

	public static final float computeSphericalDistanceToCellBorder( final int referenceCenter1, final int referenceCenter2, final int referenceCenter3, final int x, final int y, final int z, final int[][] cellBorder, final int dimX, final int dimY, final int dimZ )
	{

		float minOutsideDistance = 10000;
		int outsideIdX = 0, outsideIdY = 0, outsideIdZ = 0;

		for ( int k = 0; k < dimZ; k++ )
		{
			for ( int j = 0; j < dimY; j++ )
			{
				for ( int i = 0; i < dimX; i++ )
				{
					if ( cellBorder[ k ][ j * dimX + i ] > 0 )
					{
						final float distance = ( float ) Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) + Math.pow( z - k, 2. ) );
						if ( distance < minOutsideDistance )
						{
							minOutsideDistance = distance;
							outsideIdX = i;
							outsideIdY = j;
							outsideIdZ = k;
						}
					}
				}
			}
		}
		return ( float ) ( Math.sqrt( Math.pow( referenceCenter1 - outsideIdX, 2. ) + Math.pow( referenceCenter2 - outsideIdY, 2. ) + Math.pow( referenceCenter3 - outsideIdZ, 2. ) ) );
	}

	public static final float computeSphericalDistanceToCellBorder( final int referenceCenter1, final int referenceCenter2, final int referenceCenter3, final int x, final int y, final int z, final int[][] cellBorder, final int[][] innerBorder, final int dimX, final int dimY, final int dimZ )
	{
		float minInsideDistance = 10000, minOutsideDistance = 10000;
		int insideIdX = 0, insideIdY = 0, insideIdZ = 0, outsideIdX = 0, outsideIdY = 0, outsideIdZ = 0;

		for ( int k = 0; k < dimZ; k++ )
		{
			for ( int j = 0; j < dimY; j++ )
			{
				for ( int i = 0; i < dimX; i++ )
				{
					if ( cellBorder[ k ][ j * dimX + i ] > 0 )
					{
						final float distance = ( float ) Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) + Math.pow( z - k, 2. ) );
						if ( distance < minOutsideDistance )
						{
							minOutsideDistance = distance;
							outsideIdX = i;
							outsideIdY = j;
							outsideIdZ = k;
						}
					}
				}
			}
		}

		for ( int k = 0; k < dimZ; k++ )
		{
			for ( int j = 0; j < dimY; j++ )
			{
				for ( int i = 0; i < dimX; i++ )
				{
					if ( innerBorder[ k ][ j * dimX + i ] > 0 )
					{
						final float distance = ( float ) Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) + Math.pow( z - k, 2. ) );
						if ( distance < minInsideDistance )
						{
							minInsideDistance = distance;
							insideIdX = i;
							insideIdY = j;
							insideIdZ = k;
						}
					}
				}
			}
		}

		return ( float ) ( Math.sqrt( Math.pow( insideIdX - outsideIdX, 2. ) + Math.pow( insideIdY - outsideIdY, 2. ) + Math.pow( insideIdZ - outsideIdZ, 2. ) ) );
	}

	public static final float computeActualDistanceToCellBorder( final int x, final int y, final int z, final int[][] cellBorder, final int dimX, final int dimY, final int dimZ )
	{

		float minOutsideDistance = 10000;

		for ( int k = 0; k < dimZ; k++ )
		{
			for ( int j = 0; j < dimY; j++ )
			{
				for ( int i = 0; i < dimX; i++ )
				{
					if ( cellBorder[ k ][ j * dimX + i ] > 0 )
					{
						final float distance = ( float ) Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) + Math.pow( z - k, 2. ) );
						if ( distance < minOutsideDistance )
						{
							minOutsideDistance = distance;
						}
					}
				}
			}
		}
		return minOutsideDistance;
	}

	public static final float computeActualDistanceToCellBorder( final int x, final int y, final int z, final int[][] cellBorder, final int[][] innerBorder, final int dimX, final int dimY, final int dimZ )
	{
		float minDistance = 10000;

		for ( int k = 0; k < dimZ; k++ )
		{
			for ( int j = 0; j < dimY; j++ )
			{
				for ( int i = 0; i < dimX; i++ )
				{
					if ( cellBorder[ k ][ j * dimX + i ] > 0 )
					{
						final float distance = ( float ) Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) + Math.pow( z - k, 2. ) );
						if ( distance < minDistance )
						{
							minDistance = distance;
						}
					}
				}
			}
		}

		for ( int k = 0; k < dimZ; k++ )
		{
			for ( int j = 0; j < dimY; j++ )
			{
				for ( int i = 0; i < dimX; i++ )
				{
					if ( innerBorder[ k ][ j * dimX + i ] > 0 )
					{
						final float distance = ( float ) Math.sqrt( Math.pow( x - i, 2. ) + Math.pow( y - j, 2. ) + Math.pow( z - k, 2. ) );
						if ( distance < minDistance )
						{
							minDistance = distance;
						}
					}
				}
			}
		}

		return minDistance;
	}
}
