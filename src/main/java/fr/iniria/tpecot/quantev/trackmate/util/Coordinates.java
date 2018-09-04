package fr.iniria.tpecot.quantev.trackmate.util;

public class Coordinates
{

	/**
	 * Extract event coordinates from images.
	 *
	 * @param input
	 * @param x
	 * @param y
	 * @param z
	 * @param rad
	 * @param orientation
	 * @param segment
	 * @param considerIntensity
	 * @param radius
	 * @param angle
	 * @param depth
	 * @param intensity
	 * @param distanceToBorder
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @param i
	 */
	public static final void extractCylindricalCoordinates( final float[][] input, final int x, final int y, final int z, final float rad, final float orientation, final float segment, final boolean considerIntensity,
			final float[] radius, final float[] angle, final float[] depth, final float[] intensity, final float[] distanceToBorder,
			final int dimX, final int dimY, final int dimZ, final int i )
	{

		if ( input[ z ][ y * dimX + x ] > 0 )
		{
			angle[ i ] = orientation;
			if ( considerIntensity )
			{
				intensity[ i ] = input[ z ][ y * dimX + x ];
			}
			else
			{
				intensity[ i ] = 1;
			}
			depth[ i ] = z;
			radius[ i ] = rad;
			distanceToBorder[ i ] = segment;
		}
	}

	/**
	 * Compute angle and distance.
	 *
	 * @param angle
	 * @param distance
	 * @param intensitySum
	 * @param mask
	 * @param distanceForAngle
	 * @param xRef
	 * @param yRef
	 * @param angleUnit
	 * @param dimX
	 * @param dimY
	 * @param nbTrajectories
	 */
	public static final void computeAngleAndDistance( final float[] angle, final float[] distance, final float[] intensitySum, final int[][] mask, final float[] distanceForAngle, final int xRef, final int yRef, final double angleUnit, final int dimX, final int dimY, final int nbTrajectories )
	{
		int i = 0;
		for ( int y = 0; y < dimY; y++ )
		{
			for ( int x = 0; x < dimX; x++ )
			{
				if ( ( mask[ 0 ][ y * dimX + x ] > 0 ) && ( intensitySum[ y * dimX + x ] > 0 ) )
				{
					double theta = Math.atan2( yRef - y, x - xRef );
					if ( theta < 0 )
					{
						theta += ( 2 * Math.PI );
					}
					if ( theta > ( 2 * Math.PI ) )
					{
						theta -= ( 2 * Math.PI );
					}
					angle[ i ] = ( float ) theta;
					final CubicInterpolator interpolator = new CubicInterpolator();
					final double[] coordinates = new double[ 4 ];
					final int lowerUnit = ( int ) Math.floor( theta / angleUnit );
					int lowerUnitM = lowerUnit - 1, upperUnit = lowerUnit + 1, upperUnitP = lowerUnit + 2;
					if ( lowerUnitM < 0 )
					{
						lowerUnitM = nbTrajectories - 1;
					}
					if ( upperUnit > ( nbTrajectories - 1 ) )
					{
						upperUnit -= nbTrajectories;
					}
					if ( upperUnitP > ( nbTrajectories - 1 ) )
					{
						upperUnitP -= nbTrajectories;
					}
					coordinates[ 0 ] = distanceForAngle[ nbTrajectories + lowerUnitM ];
					coordinates[ 1 ] = distanceForAngle[ nbTrajectories + lowerUnit ];
					coordinates[ 2 ] = distanceForAngle[ nbTrajectories + upperUnit ];
					coordinates[ 3 ] = distanceForAngle[ nbTrajectories + upperUnitP ];
					final double lowerMultiplier = theta / ( float ) angleUnit - lowerUnit, xCoord = interpolator.getValue( coordinates, lowerMultiplier );
					coordinates[ 0 ] = distanceForAngle[ 2 * nbTrajectories + lowerUnitM ];
					coordinates[ 1 ] = distanceForAngle[ 2 * nbTrajectories + lowerUnit ];
					coordinates[ 2 ] = distanceForAngle[ 2 * nbTrajectories + upperUnit ];
					coordinates[ 3 ] = distanceForAngle[ 2 * nbTrajectories + upperUnitP ];
					final double yCoord = interpolator.getValue( coordinates, lowerMultiplier );
					final double currentDistance = Math.sqrt( Math.pow( xRef - xCoord, 2. ) + Math.pow( yRef - yCoord, 2. ) );
					distance[ i ] = ( float ) currentDistance;
					i++;
				}
			}
		}
	}

	/**
	 * Compute angle and distance.
	 *
	 * @param angle
	 * @param distance
	 * @param intensity
	 * @param mask
	 * @param distanceForAngle
	 * @param xRef
	 * @param yRef
	 * @param angleUnit
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @param nbTrajectories
	 */
	public static final void computeAngleAndDistance( final float[] angle, final float[] distance, final float[][] intensity, final int[][] mask, final float[] distanceForAngle, final int xRef, final int yRef, final double angleUnit, final int dimX, final int dimY, final int dimZ, final int nbTrajectories )
	{
		int i = 0;
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				for ( int x = 0; x < dimX; x++ )
				{
					if ( ( mask[ z ][ y * dimX + x ] > 0 ) && ( intensity[ z ][ y * dimX + x ] > 0 ) )
					{
						float theta = ( float ) Math.atan2( yRef - y, x - xRef );
						if ( theta < 0 )
						{
							theta += ( 2 * Math.PI );
						}
						if ( theta > ( 2 * Math.PI ) )
						{
							theta -= ( 2 * Math.PI );
						}
						angle[ i ] = theta;
						final int lowerUnit = ( int ) Math.floor( theta / angleUnit );
						final float lowerMultiplier = theta / ( float ) angleUnit - lowerUnit;
						final double currentDistance = Math.sqrt( Math.pow( xRef - ( distanceForAngle[ nbTrajectories + lowerUnit ] * lowerMultiplier + distanceForAngle[ nbTrajectories + lowerUnit + 1 ] * ( 1 - lowerMultiplier ) ), 2. )
								+ Math.pow( xRef - ( distanceForAngle[ 2 * nbTrajectories + lowerUnit ] * lowerMultiplier + distanceForAngle[ 2 * nbTrajectories + lowerUnit + 1 ] * ( 1 - lowerMultiplier ) ), 2. ) );
						distance[ i ] = ( float ) currentDistance;
						i++;
					}
				}
			}
		}
	}

	/**
	 * Extract event coordinates from images.
	 *
	 * @param input
	 * @param x
	 * @param y
	 * @param z
	 * @param rad
	 * @param orientation
	 * @param segment
	 * @param considerIntensity
	 * @param radius
	 * @param angle1
	 * @param angle2
	 * @param intensity
	 * @param distanceToBorder
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @param i
	 */
	public static final void extractSphericalCoordinates( final float[][] input, final int x, final int y, final int z, final float rad, final float orientation, final float segment, final boolean considerIntensity,
			final float[] radius, final float[] angle1, final float[] angle2, final float[] intensity, final float[] distanceToBorder,
			final int dimX, final int dimY, final int dimZ, final int i )
	{

		angle1[ i ] = orientation;
		if ( considerIntensity )
		{
			intensity[ i ] = input[ z ][ y * dimX + x ];
		}
		else
		{
			intensity[ i ] = 1;
		}
		angle2[ i ] = ( float ) Math.acos( z / Math.sqrt( rad * rad + z * z ) );
		radius[ i ] = ( float ) Math.sqrt( rad * rad + z * z );
		distanceToBorder[ i ] = segment;
	}

	/**
	 * Extract event coordinates from images.
	 *
	 * @param input
	 * @param x
	 * @param y
	 * @param z
	 * @param segment
	 * @param considerIntensity
	 * @param distance
	 * @param intensity
	 * @param unity
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @param i
	 */
	public static final void extractDistanceCoordinates( final float[][] input, final int x, final int y, final int z, final float segment, final boolean considerIntensity,
			final float[] distance, final float[] intensity, final float[] unity, final int dimX, final int dimY, final int dimZ, final int i )
	{

		distance[ i ] = segment;
		if ( considerIntensity )
		{
			intensity[ i ] = input[ z ][ y * dimX + x ];
		}
		else
		{
			intensity[ i ] = 1;
		}
		unity[ i ] = 1;
	}

	/**
	 * Extract event coordinates from images
	 *
	 * @param input
	 * @param x
	 * @param y
	 * @param z
	 * @param rad
	 * @param orientation
	 * @param segment
	 * @param considerIntensity
	 * @param xTab
	 * @param yTab
	 * @param depth
	 * @param intensity
	 * @param distanceToBorder
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @param i
	 */
	public static final void extractCartesianCoordinates( final float[][] input, final int x, final int y, final int z, final float rad, final float orientation, final float segment, final boolean considerIntensity,
			final float[] xTab, final float[] yTab, final float[] depth, final float[] intensity, final float[] distanceToBorder,
			final int dimX, final int dimY, final int dimZ, final int i )
	{

		xTab[ i ] = rad * ( float ) Math.cos( orientation );
		yTab[ i ] = rad * ( float ) Math.sin( orientation );
		if ( considerIntensity )
		{
			intensity[ i ] = input[ z ][ y * dimX + x ];
		}
		else
		{
			intensity[ i ] = 1;
		}
		depth[ i ] = z;
		distanceToBorder[ i ] = segment;
	}
}
