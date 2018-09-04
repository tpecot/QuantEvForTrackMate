package fr.iniria.tpecot.quantev.trackmate.util;

/**
 * Morphological mathematics
 */
public class Morphology
{

	/**
	 * Structuring element
	 *
	 * @param radiusXY
	 * @param radiusZ
	 * @return
	 */
	public static final int[][] makeStrElt( final int radiusXY, final int radiusZ )
	{
		final int StrEltDimXY = 2 * radiusXY + 1, StrEltDimZ = 2 * radiusZ + 1;
		final int[][] StrElt = new int[ StrEltDimZ ][ StrEltDimXY * StrEltDimXY ];
		for ( int dz = -radiusZ; dz <= radiusZ; dz++ )
		{
			for ( int dy = -radiusXY; dy <= radiusXY; dy++ )
			{
				for ( int dx = -radiusXY; dx <= radiusXY; dx++ )
				{
					if ( Math.sqrt( Math.pow( dx, 2d ) + Math.pow( dy, 2d ) + Math.pow( dz, 2d ) ) < radiusXY )
					{
						StrElt[ dz + radiusZ ][ ( dy + radiusXY ) * StrEltDimXY + dx + radiusXY ] = 1;
					}
				}
			}
		}
		return StrElt;
	}

	/**
	 * Dilation
	 *
	 * @param data
	 * @param StrElt
	 * @param radiusXY
	 * @param radiusZ
	 * @param dimX
	 * @param dimY
	 * @param dimZ
	 * @return
	 */
	public static final int[][] dilate( final int[][] data, final int[][] StrElt, final int radiusXY, final int radiusZ, final int dimX, final int dimY, final int dimZ )
	{
		final int[][] result = new int[ dimZ ][ dimX * dimY ];
		final int StrEltDimXY = 2 * radiusXY + 1, StrEltDimZ = 2 * radiusZ + 1;
		for ( int z = 0; z < dimZ; z++ )
		{
			for ( int y = 0; y < dimY; y++ )
			{
				for ( int x = 0; x < dimX; x++ )
				{
					if ( data[ z ][ y * dimX + x ] > 0 )
					{
						for ( int w = 0; w < StrEltDimZ; w++ )
						{
							for ( int v = 0; v < StrEltDimXY; v++ )
							{
								for ( int u = 0; u < StrEltDimXY; u++ )
								{
									if ( ( ( x + u - radiusXY ) >= 0 ) && ( ( x + u - radiusXY ) < dimX ) && ( ( y + v - radiusXY ) >= 0 ) && ( ( y + v - radiusXY ) < dimY ) && ( ( z + w - radiusZ ) >= 0 ) && ( ( z + w - radiusZ ) < dimZ ) )
									{
										if ( StrElt[ w ][ v * StrEltDimXY + u ] > 0.0001 )
										{
											result[ z + w - radiusZ ][ ( y + v - radiusXY ) * dimX + x + u - radiusXY ] = 1;
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
}
