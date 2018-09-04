package fr.iniria.tpecot.quantev.trackmate.util;

public class BesselFunctions
{

	/**
	 * Modified Bessel function.
	 *
	 * @param N
	 * @param X
	 * @return
	 */
	public static final double bessel_i( final int N, final double X )
	{
		final double BIGNO = 1e10, BIGNI = 1e-10;
		final int IACC = 40, M = ( int ) ( 2 * ( ( N + Math.floor( Math.sqrt( IACC * N ) ) ) ) );
		if ( N == 0 )
			return ( bessel_i0( X ) );
		if ( N == 1 )
			return ( bessel_i1( X ) );
		if ( X == 0.0 )
			return 0.0;
		final double TOX = 2.0 / X;
		double BIP = 0.0, BI = 1.0, BSI = 0.0;
		for ( int J = M; J > 0; J-- )
		{
			final double BIM = BIP + J * TOX * BI;
			BIP = BI;
			BI = BIM;
			if ( Math.abs( BI ) > BIGNO )
			{
				BI = BI * BIGNI;
				BIP = BIP * BIGNI;
				BSI = BSI * BIGNI;
			}
			if ( J == N )
				BSI = BIP;
		}
		return ( BSI * bessel_i0( X ) / BI );
	}

	/**
	 * Modified Bessel function of order 0.
	 *
	 * @param X
	 * @return
	 */
	public static final double bessel_i0( final double X )
	{
		final double P1 = 1.0, P2 = 3.5156229, P3 = 3.0899424, P4 = 1.2067492,
				P5 = 0.2659732, P6 = 0.360768e-1, P7 = 0.45813e-2,
				Q1 = 0.39894228, Q2 = 0.1328592e-1, Q3 = 0.225319e-2,
				Q4 = -0.157565e-2, Q5 = 0.916281e-2, Q6 = -0.2057706e-1,
				Q7 = 0.2635537e-1, Q8 = -0.1647633e-1, Q9 = 0.392377e-2,
				AX = Math.abs( X );
		if ( AX < 3.75 )
		{
			final double Y = ( X / 3.75 ) * ( X / 3.75 );
			return ( P1 + Y * ( P2 + Y * ( P3 + Y * ( P4 + Y * ( P5 + Y * ( P6 + Y * P7 ) ) ) ) ) );
		}
		else
		{
			final double Y = 3.75 / AX, BX = Math.exp( AX ) / Math.sqrt( AX );
			return ( BX * ( Q1 + Y * ( Q2 + Y * ( Q3 + Y * ( Q4 + Y * ( Q5 + Y * ( Q6 + Y * ( Q7 + Y * ( Q8 + Y * Q9 ) ) ) ) ) ) ) ) );
		}
	}

	/**
	 * Modified Bessel function of order 1.
	 *
	 * @param X
	 * @return
	 */
	private static final double bessel_i1( final double X )
	{
		final double P1 = 0.5, P2 = 0.87890594, P3 = 0.51498869, P4 = 0.15084934,
				P5 = 0.2658733e-1, P6 = 0.301532e-2, P7 = 0.32411e-3,
				Q1 = 0.39894228, Q2 = -0.3988024e-1, Q3 = -0.362018e-2,
				Q4 = 0.163801e-2, Q5 = -0.1031555e-1, Q6 = 0.2282967e-1,
				Q7 = -0.2895312e-1, Q8 = 0.1787654e-1, Q9 = -0.420059e-2,
				AX = Math.abs( X );
		if ( AX < 3.75 )
		{
			final double Y = ( X / 3.75 ) * ( X / 3.75 );
			return ( X * ( P1 + Y * ( P2 + Y * ( P3 + Y * ( P4 + Y * ( P5 + Y * ( P6 + Y * P7 ) ) ) ) ) ) );
		}
		else
		{
			final double Y = 3.75 / AX, BX = Math.exp( AX ) / Math.sqrt( AX );
			return ( BX * ( Q1 + Y * ( Q2 + Y * ( Q3 + Y * ( Q4 + Y * ( Q5 + Y * ( Q6 + Y * ( Q7 + Y * ( Q8 + Y * Q9 ) ) ) ) ) ) ) ) );
		}
	}
}
