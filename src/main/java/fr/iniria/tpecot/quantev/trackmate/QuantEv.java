package fr.iniria.tpecot.quantev.trackmate;

import java.io.File;

import static fr.iniria.tpecot.quantev.trackmate.util.BesselFunctions.*;
import jxl.Cell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;

import jsc.datastructures.*;
import jsc.tests.*;
import jsc.onesample.*;

public class QuantEv
{

	// compute bandwidth according to the rule of thumb
	float compute_rule_of_thumb_bandwidth( final float[] input )
	{

		float mu = 0;

		// Compute mean and standard deviation
		float mean = 0, std = 0;
		for ( int i = 0; i < input.length; i++ )
		{
			mean += input[ i ];
		}
		mean /= input.length - 1;
		for ( int i = 0; i < input.length; i++ )
		{
			std += ( float ) ( Math.pow( ( ( double ) input[ i ] - mean ), 2. ) );
		}
		std = ( float ) Math.sqrt( ( double ) std / ( double ) ( input.length - 1 ) );
		mu = ( float ) ( Math.pow( ( 4.0 * Math.pow( std, 5. ) ) / ( 3.0 * input.length ), 1.0 / 5.0 ) );

		return mu;
	}

	// compute histogram for radius
	float[] computeRadiusHistogram( final float[] input, final float[] weight, final float[] length, final int nbBins )
	{

		final float[] histogram = new float[ nbBins ];
		final int[] nbEltsPerBin = new int[ nbBins ];
		float maxValue = 0, lengthMean = 0, lengthVar = 0, nbEvtsPerExp = 0;
		boolean takeLengthIntoAccount = false;
		for ( int i = 0; i < input.length; i++ )
		{
			lengthMean += length[ i ];
			nbEvtsPerExp += 1.;
		}
		lengthMean /= nbEvtsPerExp;
		for ( int i = 0; i < input.length; i++ )
		{
			lengthVar += Math.pow( length[ i ] - lengthMean, 2. );
		}
		lengthVar /= nbEvtsPerExp;
		if ( lengthVar > 1. )
		{
			takeLengthIntoAccount = true;
		}

		if ( !takeLengthIntoAccount )
		{
			for ( int i = 0; i < input.length; i++ )
			{
				if ( input[ i ] > maxValue )
				{
					maxValue = input[ i ];
				}
			}
		}
		for ( int i = 0; i < input.length; i++ )
		{
			int currentBin;
			if ( takeLengthIntoAccount )
			{
				if ( ( !Float.isNaN( input[ i ] ) ) && ( !Float.isNaN( length[ i ] ) ) )
				{
					currentBin = Math.round( ( nbBins - 1 ) * input[ i ] / length[ i ] );
					if ( !Float.isNaN( currentBin ) )
					{
						if ( currentBin == nbBins )
						{
							currentBin = 0;
						}
						if ( !Float.isNaN( weight[ i ] ) && !Float.isNaN( length[ i ] ) )
						{
							histogram[ currentBin ] += weight[ i ];
							nbEltsPerBin[ currentBin ]++;
						}
					}
				}
			}
			else
			{
				if ( !Float.isNaN( input[ i ] ) )
				{
					currentBin = Math.round( ( nbBins - 1 ) * input[ i ] / maxValue );
					if ( !Float.isNaN( currentBin ) )
					{
						if ( currentBin == nbBins )
						{
							currentBin = 0;
						}
						if ( !Float.isNaN( weight[ i ] ) && !Float.isNaN( length[ i ] ) && ( weight[ i ] < 1000000 ) )
						{
							histogram[ currentBin ] += weight[ i ];
							nbEltsPerBin[ currentBin ]++;
						}
					}
				}
			}
		}

		for ( int i = 0; i < nbBins; i++ )
		{
			if ( nbEltsPerBin[ i ] > 0 )
			{
				histogram[ i ] /= nbEltsPerBin[ i ];
			}
		}
		return histogram;
	}

	// compute histogram for radius
	float[] computeDepthHistogram( final float[] input, final float[] weight, final float[] length, final int nbBins, final float maxValue, float[] interpolatedHistogram )
	{

		int actualNbBins;
		if ( nbBins > ( maxValue + 1. ) )
		{
			actualNbBins = ( int ) Math.ceil( maxValue + 1 );
		}
		else
		{
			actualNbBins = nbBins;
		}

		final float[] histogram = new float[ actualNbBins ];
		float[] outputHistogram = new float[ nbBins ];
		final int[] nbEltsPerBin = new int[ nbBins ];
		for ( int i = 0; i < input.length; i++ )
		{
			if ( !Float.isNaN( input[ i ] ) )
			{
				int currentBin = Math.round( ( actualNbBins - 1 ) * input[ i ] / maxValue );
				if ( currentBin == actualNbBins )
				{
					currentBin = 0;
				}
				if ( !Float.isNaN( weight[ i ] ) && !Float.isNaN( length[ i ] ) && ( weight[ i ] < 1000000 ) )
				{
					histogram[ currentBin ] += weight[ i ];
					nbEltsPerBin[ currentBin ]++;
				}
			}
		}
		final double totalSum = 0.;
		for ( int i = 0; i < actualNbBins; i++ )
		{
			if ( nbEltsPerBin[ i ] > 0 )
			{
				histogram[ i ] /= nbEltsPerBin[ i ];
			}
		}

		if ( actualNbBins != nbBins )
		{
			final int translationBin = ( int ) ( nbBins / ( 2 * ( float ) actualNbBins ) );
			for ( int i = 0; i < translationBin; i++ )
			{
				interpolatedHistogram[ i ] = histogram[ 0 ];
			}
			for ( int i = 0; i < nbBins; i++ )
			{
				outputHistogram[ i ] = histogram[ ( int ) ( ( float ) i * ( float ) actualNbBins / nbBins ) ];
				final int lowerBound = ( int ) Math.floor( ( float ) i * ( float ) actualNbBins / nbBins ),
						higherBound = ( int ) Math.ceil( ( float ) i * ( float ) actualNbBins / nbBins );
				if ( i < ( nbBins - translationBin ) )
				{
					if ( lowerBound == higherBound )
					{
						interpolatedHistogram[ i + translationBin ] = histogram[ lowerBound ];
					}
					else
					{
						if ( higherBound < actualNbBins )
						{
							interpolatedHistogram[ i + translationBin ] = ( 1 - Math.abs( i - ( float ) lowerBound * ( float ) nbBins / actualNbBins ) / ( ( float ) nbBins / ( float ) actualNbBins ) ) * ( histogram[ lowerBound ] )
									+ ( 1 - Math.abs( i - ( float ) higherBound * ( float ) nbBins / actualNbBins ) / ( ( float ) nbBins / ( float ) actualNbBins ) ) * ( histogram[ higherBound ] );
						}
						else
						{
							interpolatedHistogram[ i + translationBin ] = histogram[ actualNbBins - 1 ];
						}
					}
				}
				else
				{
					interpolatedHistogram[ i ] = histogram[ actualNbBins - 1 ];
				}
			}
		}
		else
		{
			outputHistogram = histogram;
			interpolatedHistogram = histogram;
		}

		return outputHistogram;
	}

	// compute Gaussian kernel for density estimation
	float[] computeGaussianKernel( final int nbBins, final float bandwidth, final float maxValue )
	{

		final float[] GaussianKernel = new float[ nbBins ];
		for ( int i = 0; i < nbBins; i++ )
		{
			final float correspondingValue = ( ( ( i ) - ( float ) ( nbBins ) / 2 ) * maxValue / ( nbBins ) );
			GaussianKernel[ i ] = ( float ) ( Math.exp( -Math.pow( correspondingValue / bandwidth, 2. ) / 2. ) / ( bandwidth * Math.sqrt( 2 * Math.PI ) ) );
		}

		return GaussianKernel;
	}

	// compute density
	float[] computeDensity( final int nbBins, final float[] histogram, final float[] GaussianKernel )
	{
		final float[] density = new float[ nbBins ];
		for ( int i = 0; i < nbBins; i++ )
		{
			float GaussianNormalization = 0;
			for ( int j = 0; j < nbBins; j++ )
			{
				final int currentIndex = i + j - nbBins / 2;
				if ( ( currentIndex >= 0 ) && ( currentIndex < nbBins ) )
				{
					density[ i ] += histogram[ currentIndex ] * GaussianKernel[ j ];
					GaussianNormalization += GaussianKernel[ j ];
				}
			}
			density[ i ] /= GaussianNormalization;
		}

		return density;
	}

	// compute histogram for radius
	float[] computeSymmetricalHistogram( final float[] input, final float[] weight, final float[] length, final int nbBins, final float minValue, final float maxValue, float[] interpolatedHistogram )
	{

		int actualNbBins;
		if ( nbBins > ( 2 * maxValue + 1. ) )
		{
			actualNbBins = ( int ) Math.ceil( 2 * maxValue + 1 );
		}
		else
		{
			actualNbBins = nbBins;
		}

		final float[] histogram = new float[ actualNbBins ];
		float[] outputHistogram = new float[ nbBins ];
		final int[] nbEltsPerBin = new int[ nbBins ];
		for ( int i = 0; i < input.length; i++ )
		{
			if ( !Float.isNaN( input[ i ] ) )
			{
				int currentBin = Math.round( ( actualNbBins - 1 ) * ( input[ i ] - minValue ) / ( 2 * maxValue ) );
				if ( currentBin == actualNbBins )
				{
					currentBin = 0;
				}
				if ( !Float.isNaN( weight[ i ] ) && !Float.isNaN( length[ i ] ) && ( weight[ i ] < 1000000 ) )
				{
					histogram[ currentBin ] += weight[ i ];
					nbEltsPerBin[ i ]++;
				}
			}
		}

		final double totalSum = 0.;
		for ( int i = 0; i < actualNbBins; i++ )
		{
			if ( nbEltsPerBin[ i ] > 0 )
			{
				histogram[ i ] /= nbEltsPerBin[ i ];
			}
		}

		if ( actualNbBins != nbBins )
		{
			final int translationBin = ( int ) ( nbBins / ( 2 * ( float ) actualNbBins ) );
			for ( int i = 0; i < translationBin; i++ )
			{
				interpolatedHistogram[ i ] = histogram[ 0 ];
			}
			for ( int i = 0; i < nbBins; i++ )
			{
				outputHistogram[ i ] = histogram[ ( int ) ( ( float ) i * ( float ) actualNbBins / nbBins ) ];
				final int lowerBound = ( int ) Math.floor( ( float ) i * ( float ) actualNbBins / nbBins ),
						higherBound = ( int ) Math.ceil( ( float ) i * ( float ) actualNbBins / nbBins );
				if ( i < ( nbBins - translationBin ) )
				{
					if ( lowerBound == higherBound )
					{
						interpolatedHistogram[ i + translationBin ] = histogram[ lowerBound ];
					}
					else
					{
						if ( higherBound < actualNbBins )
						{
							interpolatedHistogram[ i + translationBin ] = ( 1 - Math.abs( i - ( float ) lowerBound * ( float ) nbBins / actualNbBins ) / ( ( float ) nbBins / ( float ) actualNbBins ) ) * ( histogram[ lowerBound ] )
									+ ( 1 - Math.abs( i - ( float ) higherBound * ( float ) nbBins / actualNbBins ) / ( ( float ) nbBins / ( float ) actualNbBins ) ) * ( histogram[ higherBound ] );
						}
						else
						{
							interpolatedHistogram[ i + translationBin ] = histogram[ actualNbBins - 1 ];
						}
					}
				}
				else
				{
					interpolatedHistogram[ i ] = histogram[ actualNbBins - 1 ];
				}
			}
		}
		else
		{
			outputHistogram = histogram;
			interpolatedHistogram = histogram;
		}

		return outputHistogram;
	}

	// compute Gaussian kernel for density estimation
	float[] computeSymmetricalGaussianKernel( final int nbBins, final float bandwidth, final float maxValue )
	{

		final float[] GaussianKernel = new float[ nbBins ];
		for ( int i = 0; i < nbBins; i++ )
		{
			final float correspondingValue = ( ( ( i ) - ( float ) ( nbBins ) / 2 ) * ( 2 * maxValue + 1 ) / ( nbBins ) );
			GaussianKernel[ i ] = ( float ) ( Math.exp( -Math.pow( correspondingValue / bandwidth, 2. ) / 2. ) / ( bandwidth * Math.sqrt( 2 * Math.PI ) ) );
		}

		return GaussianKernel;
	}

	// ! Rule of Thumb estimation of the concentration for von Mises
	// distribution kernel
	/**
	 * \return the concentration for the van mises distribution \param S a
	 * vector of angles
	 **/
	float compute_rule_of_thumb_concentration( final float[] input )
	{

		float mu;

		// Compute R
		double kappaRef = 0.;
		for ( int p = 1; p <= 3; p++ )
		{
			double sc = 0, ss = 0;
			final double n = input.length;
			double R = 0.;
			for ( int i = 0; i < input.length; i++ )
			{
				final double theta = input[ i ];
				sc += Math.cos( p * theta );
				ss += Math.sin( p * theta );
			}
			final double m = Math.atan2( ss, sc );
			for ( int i = 0; i < input.length; i++ )
			{
				final double theta = input[ i ];
				R += Math.cos( p * theta - m );
			}
			R /= n;
			// Estimate the concentration kappa of the data using a newton
			// method
			// http://en.wikipedia.org/wiki/Von_Mises%E2%80%93Fisher_distribution
			double kappa = R * ( 2.0 - R * R ) / ( 1.0 - R * R ), dkappa = 1.0;
			for ( int i = 0; i < 1000 && Math.abs( dkappa ) > 1e-12; i++ )
			{
				final double A = bessel_i( p, kappa ) / bessel_i0( kappa );
				dkappa = ( A - R ) / ( 1.0 - A * A - A / kappa );
				kappa -= dkappa;
			}
			if ( kappa > kappaRef )
			{
				kappaRef = kappa;
			}
		}
		// from (Taylor, 2008) deduce the plugin concentration for smoothing
		mu = ( float ) Math.pow( ( 3.0 * ( input.length ) * kappaRef * kappaRef * bessel_i( 2, 2.0 * kappaRef ) ) / ( 4.0 * Math.sqrt( Math.PI ) * bessel_i0( kappaRef ) * bessel_i0( kappaRef ) ), 2.0 / 5.0 );
		if ( Float.isNaN( mu ) )
		{
			mu = 0;
		}

		return mu;

	}

	// compute histogram for theta distribution over several experiments
	float[] computeCircularHistogram( final float[] angle, final float[] weight, final float[] length, final int nbBins )
	{

		final float[] angleHistogram = new float[ nbBins ];
		final int[] nbEltsPerBin = new int[ nbBins ];
		for ( int i = 0; i < angle.length; i++ )
		{
			if ( !Float.isNaN( angle[ i ] ) )
			{
				int currentBin = ( int ) Math.round( nbBins * angle[ i ] / ( 2 * Math.PI ) );
				if ( currentBin == nbBins )
				{
					currentBin = 0;
				}
				if ( !Float.isNaN( weight[ i ] ) && !Float.isNaN( length[ i ] ) )
				{
					angleHistogram[ currentBin ] += weight[ i ];
					nbEltsPerBin[ currentBin ]++;
				}
			}
		}

		for ( int i = 0; i < nbBins; i++ )
		{
			if ( nbEltsPerBin[ i ] > 0 )
			{
				angleHistogram[ i ] /= nbEltsPerBin[ i ];
			}
		}
		return angleHistogram;
	}

	// compute von Mises kernel for circular density estimation
	float[] computeVonMisesKernel( final int nbBins, float concentration )
	{

		final float[] vonMisesKernel = new float[ nbBins ];
		if ( concentration > 709. )
		{
			concentration = 709;
		}
		for ( int i = 0; i < nbBins; i++ )
		{
			final float correspondingValue = ( ( i - ( float ) nbBins / 2 ) * 2 * ( float ) Math.PI / nbBins );
			vonMisesKernel[ i ] = ( float ) ( Math.exp( concentration * Math.cos( correspondingValue ) ) / ( 2 * Math.PI * bessel_i0( concentration ) ) );
		}

		return vonMisesKernel;
	}

	// compute circular density
	float[] computeCircularDensity( final int nbBins, final float[] histogram, final float[] vonMisesKernel )
	{

		final float[] density = new float[ nbBins ];
		float totalDensitySum = 0;
		for ( int i = 0; i < nbBins; i++ )
		{
			for ( int j = 0; j < nbBins; j++ )
			{
				int currentIndex = i + j - nbBins / 2;
				if ( currentIndex < 0 )
				{
					currentIndex += nbBins;
				}
				if ( currentIndex >= nbBins )
				{
					currentIndex -= nbBins;
				}
				density[ i ] += histogram[ currentIndex ] * vonMisesKernel[ j ];
			}
			totalDensitySum += density[ i ];
		}
		for ( int i = 0; i < nbBins; i++ )
		{
			density[ i ] = density[ i ] / totalDensitySum;
		}

		return density;
	}

	// export histograms and densities as xls file
	void exportCylindricalFiles( final float[] histogram1, final float[] histogram2, final float[] histogram3, final float[] density1, final float[] density2, final float[] density3, final File f, final Sequence seq )
	{
		if ( f.exists() )
		{
			WritableWorkbook wb = null;
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet depthHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			WritableSheet depthDensitySheet = null;
			int radiusRows = 0, angleRows = 0, depthRows = 0, currentColumn = 0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite( f );
				if ( wb.getNumberOfSheets() != 6 )
				{
					MessageDialog.showDialog( "The excel output file should contain 6 worksheets: Radius histogram, Polar angle histogram, Depth histogram, Radius density, Polar angle density and Depth density." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				radiusHistogramSheet = wb.getSheet( 0 );
				angleHistogramSheet = wb.getSheet( 1 );
				depthHistogramSheet = wb.getSheet( 2 );
				radiusDensitySheet = wb.getSheet( 3 );
				angleDensitySheet = wb.getSheet( 4 );
				depthDensitySheet = wb.getSheet( 5 );
				currentColumn = radiusHistogramSheet.getColumns();
				if ( ( angleHistogramSheet.getColumns() != currentColumn ) || ( depthHistogramSheet.getColumns() != currentColumn ) || ( radiusDensitySheet.getColumns() != currentColumn ) || ( angleDensitySheet.getColumns() != currentColumn ) || ( depthDensitySheet.getColumns() != currentColumn ) )
				{
					MessageDialog.showDialog( "The excel output file should have the same number of columns in each worksheet." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				radiusRows = radiusHistogramSheet.getRows();
				if ( radiusDensitySheet.getRows() != radiusRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Radius histogram and Radius density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				angleRows = angleHistogramSheet.getRows();
				if ( angleDensitySheet.getRows() != angleRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Polar angle histogram and Polar angle density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				depthRows = depthHistogramSheet.getRows();
				if ( depthDensitySheet.getRows() != depthRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Depth histogram and Depth density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}

			if ( histogram1.length != ( radiusRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( radiusHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( radiusDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( radiusHistogramSheet, currentColumn, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( radiusDensitySheet, currentColumn, bin + 1, density1[ bin ] );
			}
			if ( histogram2.length != ( angleRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Polar angle histogram and Polar angle density do not have the same number of rows than the number of bins for polar angle." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( angleHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( angleDensitySheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellNumber( angleHistogramSheet, currentColumn, 1, histogram2[ 0 ] );
			XLSUtil.setCellNumber( angleDensitySheet, currentColumn, 1, density2[ 0 ] );
			for ( int bin = 1; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( angleHistogramSheet, currentColumn, bin + 1, histogram2[ histogram2.length - bin ] );
				XLSUtil.setCellNumber( angleDensitySheet, currentColumn, bin + 1, density2[ histogram2.length - bin ] );
			}
			if ( histogram3.length != ( depthRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Depth histogram and Depth density do not have the same number of rows than the number of bins for depth." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( depthHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( depthDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram3.length; bin++ )
			{
				XLSUtil.setCellNumber( depthHistogramSheet, currentColumn, bin + 1, histogram3[ bin ] );
				XLSUtil.setCellNumber( depthDensitySheet, currentColumn, bin + 1, density3[ bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
		else
		{
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet depthHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			WritableSheet depthDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook( f );
				radiusHistogramSheet = XLSUtil.createNewPage( wb, "Radius histogram" );
				angleHistogramSheet = XLSUtil.createNewPage( wb, "Polar angle histogram" );
				depthHistogramSheet = XLSUtil.createNewPage( wb, "Depth histogram" );
				radiusDensitySheet = XLSUtil.createNewPage( wb, "Radius density" );
				angleDensitySheet = XLSUtil.createNewPage( wb, "Polar angle density" );
				depthDensitySheet = XLSUtil.createNewPage( wb, "Depth density" );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
			XLSUtil.setCellString( radiusHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( radiusDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( radiusHistogramSheet, 0, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( radiusDensitySheet, 0, bin + 1, density1[ bin ] );
			}
			XLSUtil.setCellString( angleHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( angleDensitySheet, 0, 0, seq.getName() );
			XLSUtil.setCellNumber( angleHistogramSheet, 0, 1, histogram2[ 0 ] );
			XLSUtil.setCellNumber( angleDensitySheet, 0, 1, density2[ 0 ] );
			for ( int bin = 1; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( angleHistogramSheet, 0, bin + 1, histogram2[ histogram2.length - bin ] );
				XLSUtil.setCellNumber( angleDensitySheet, 0, bin + 1, density2[ histogram2.length - bin ] );
			}
			XLSUtil.setCellString( depthHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( depthDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram3.length; bin++ )
			{
				XLSUtil.setCellNumber( depthHistogramSheet, 0, bin + 1, histogram3[ bin ] );
				XLSUtil.setCellNumber( depthDensitySheet, 0, bin + 1, density3[ bin ] );
			}
			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
	}

	// export histogram or density as xls file
	void exportCylindricalFiles( final float[] histogram1, final float[] histogram2, final float[] density1, final float[] density2, final File f, final Sequence seq )
	{
		if ( f.exists() )
		{
			WritableWorkbook wb = null;
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			int radiusRows = 0, angleRows = 0, currentColumn = 0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite( f );
				if ( wb.getNumberOfSheets() != 4 )
				{
					MessageDialog.showDialog( "The excel output file should contain 4 worksheets: Radius histogram, Polar angle histogram, Radius density and Polar angle density." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				radiusHistogramSheet = wb.getSheet( 0 );
				angleHistogramSheet = wb.getSheet( 1 );
				radiusDensitySheet = wb.getSheet( 2 );
				angleDensitySheet = wb.getSheet( 3 );
				currentColumn = radiusHistogramSheet.getColumns();
				if ( ( angleHistogramSheet.getColumns() != currentColumn ) || ( radiusDensitySheet.getColumns() != currentColumn ) || ( angleDensitySheet.getColumns() != currentColumn ) )
				{
					MessageDialog.showDialog( "The excel output file should have the same number of columns in each worksheet." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				radiusRows = radiusHistogramSheet.getRows();
				if ( radiusDensitySheet.getRows() != radiusRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Radius histogram and Radius density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				angleRows = angleHistogramSheet.getRows();
				if ( angleDensitySheet.getRows() != angleRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Polar angle histogram and Polar angle density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}

			if ( histogram1.length != ( radiusRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( radiusHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( radiusDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( radiusHistogramSheet, currentColumn, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( radiusDensitySheet, currentColumn, bin + 1, density1[ bin ] );
			}
			if ( histogram2.length != ( angleRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Polar angle histogram and Polar angle density do not have the same number of rows than the number of bins for polar angle." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( angleHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( angleDensitySheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellNumber( angleHistogramSheet, currentColumn, 1, histogram2[ 0 ] );
			XLSUtil.setCellNumber( angleDensitySheet, currentColumn, 1, density2[ 0 ] );
			for ( int bin = 1; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( angleHistogramSheet, currentColumn, bin + 1, histogram2[ histogram2.length - bin ] );
				XLSUtil.setCellNumber( angleDensitySheet, currentColumn, bin + 1, density2[ histogram2.length - bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
		else
		{
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angleHistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angleDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook( f );
				radiusHistogramSheet = XLSUtil.createNewPage( wb, "Radius histogram" );
				angleHistogramSheet = XLSUtil.createNewPage( wb, "Polar angle histogram" );
				radiusDensitySheet = XLSUtil.createNewPage( wb, "Radius density" );
				angleDensitySheet = XLSUtil.createNewPage( wb, "Polar angle density" );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
			XLSUtil.setCellString( radiusHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( radiusDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( radiusHistogramSheet, 0, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( radiusDensitySheet, 0, bin + 1, density1[ bin ] );
			}
			XLSUtil.setCellString( angleHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( angleDensitySheet, 0, 0, seq.getName() );
			XLSUtil.setCellNumber( angleHistogramSheet, 0, 1, histogram2[ 0 ] );
			XLSUtil.setCellNumber( angleDensitySheet, 0, 1, density2[ 0 ] );
			for ( int bin = 1; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( angleHistogramSheet, 0, bin + 1, histogram2[ histogram2.length - bin ] );
				XLSUtil.setCellNumber( angleDensitySheet, 0, bin + 1, density2[ histogram2.length - bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
	}

	// export histograms and densities as xls file
	void exportCartesianFiles( final float[] histogram1, final float[] histogram2, final float[] histogram3, final float[] density1, final float[] density2, final float[] density3, final File f, final Sequence seq )
	{
		if ( f.exists() )
		{
			WritableWorkbook wb = null;
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet zHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			WritableSheet zDensitySheet = null;
			int xRows = 0, yRows = 0, zRows = 0, currentColumn = 0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite( f );
				if ( wb.getNumberOfSheets() != 6 )
				{
					MessageDialog.showDialog( "The excel output file should contain 6 worksheets: X histogram, Y histogram, Z histogram, X density, Y density and Z density." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				xHistogramSheet = wb.getSheet( 0 );
				yHistogramSheet = wb.getSheet( 1 );
				zHistogramSheet = wb.getSheet( 2 );
				xDensitySheet = wb.getSheet( 3 );
				yDensitySheet = wb.getSheet( 4 );
				zDensitySheet = wb.getSheet( 5 );
				currentColumn = xHistogramSheet.getColumns();
				if ( ( yHistogramSheet.getColumns() != currentColumn ) || ( zHistogramSheet.getColumns() != currentColumn ) || ( xDensitySheet.getColumns() != currentColumn ) || ( yDensitySheet.getColumns() != currentColumn ) || ( zDensitySheet.getColumns() != currentColumn ) )
				{
					MessageDialog.showDialog( "The excel output file should have the same number of columns in each worksheet." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				xRows = xHistogramSheet.getRows();
				if ( xDensitySheet.getRows() != xRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets X histogram and X density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				yRows = yHistogramSheet.getRows();
				if ( yDensitySheet.getRows() != yRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Y histogram and Y density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				zRows = zHistogramSheet.getRows();
				if ( zDensitySheet.getRows() != zRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Z histogram and Z density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}

			if ( histogram1.length != ( xRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets X histogram and X density do not have the same number of rows than the number of bins for abscissa." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( xHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( xDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( xHistogramSheet, currentColumn, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( xDensitySheet, currentColumn, bin + 1, density1[ bin ] );
			}
			if ( histogram2.length != ( yRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Y histogram and Y density do not have the same number of rows than the number of bins for ordinate." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( yHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( yDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( yHistogramSheet, currentColumn, bin + 1, histogram2[ bin ] );
				XLSUtil.setCellNumber( yDensitySheet, currentColumn, bin + 1, density2[ bin ] );
			}
			if ( histogram3.length != ( zRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Z histogram and Z density do not have the same number of rows than the number of bins for height." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( zHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( zDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram3.length; bin++ )
			{
				XLSUtil.setCellNumber( zHistogramSheet, currentColumn, bin + 1, histogram3[ bin ] );
				XLSUtil.setCellNumber( zDensitySheet, currentColumn, bin + 1, density3[ bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
		else
		{
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet zHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			WritableSheet zDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook( f );
				xHistogramSheet = XLSUtil.createNewPage( wb, "X histogram" );
				yHistogramSheet = XLSUtil.createNewPage( wb, "Y histogram" );
				zHistogramSheet = XLSUtil.createNewPage( wb, "Z histogram" );
				xDensitySheet = XLSUtil.createNewPage( wb, "X density" );
				yDensitySheet = XLSUtil.createNewPage( wb, "Y density" );
				zDensitySheet = XLSUtil.createNewPage( wb, "Z density" );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
			XLSUtil.setCellString( xHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( xDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( xHistogramSheet, 0, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( xDensitySheet, 0, bin + 1, density1[ bin ] );
			}
			XLSUtil.setCellString( yHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( yDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( yHistogramSheet, 0, bin + 1, histogram2[ bin ] );
				XLSUtil.setCellNumber( yDensitySheet, 0, bin + 1, density2[ bin ] );
			}
			XLSUtil.setCellString( zHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( zDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram3.length; bin++ )
			{
				XLSUtil.setCellNumber( zHistogramSheet, 0, bin + 1, histogram3[ bin ] );
				XLSUtil.setCellNumber( zDensitySheet, 0, bin + 1, density3[ bin ] );
			}
			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
	}

	// export histogram or density as xls file
	void exportCartesianFiles( final float[] histogram1, final float[] histogram2, final float[] density1, final float[] density2, final File f, final Sequence seq )
	{
		if ( f.exists() )
		{
			WritableWorkbook wb = null;
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			int xRows = 0, yRows = 0, currentColumn = 0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite( f );
				if ( wb.getNumberOfSheets() != 4 )
				{
					MessageDialog.showDialog( "The excel output file should contain 4 worksheets: X histogram, Y histogram, X density and Y density." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				xHistogramSheet = wb.getSheet( 0 );
				yHistogramSheet = wb.getSheet( 1 );
				xDensitySheet = wb.getSheet( 2 );
				yDensitySheet = wb.getSheet( 3 );
				currentColumn = xHistogramSheet.getColumns();
				if ( ( yHistogramSheet.getColumns() != currentColumn ) || ( xDensitySheet.getColumns() != currentColumn ) || ( yDensitySheet.getColumns() != currentColumn ) )
				{
					MessageDialog.showDialog( "The excel output file should have the same number of columns in each worksheet." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				xRows = xHistogramSheet.getRows();
				if ( xDensitySheet.getRows() != xRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets X histogram and X density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				yRows = yHistogramSheet.getRows();
				if ( yDensitySheet.getRows() != yRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Y histogram and Y density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}

			if ( histogram1.length != ( xRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets X histogram and X density do not have the same number of rows than the number of bins for abscissa." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( xHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( xDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( xHistogramSheet, currentColumn, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( xDensitySheet, currentColumn, bin + 1, density1[ bin ] );
			}
			if ( histogram2.length != ( yRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Y histogram and Y density do not have the same number of rows than the number of bins for ordinate." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( yHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( yDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( yHistogramSheet, currentColumn, bin + 1, histogram2[ bin ] );
				XLSUtil.setCellNumber( yDensitySheet, currentColumn, bin + 1, density2[ bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
		else
		{
			WritableSheet xHistogramSheet = null;
			WritableSheet yHistogramSheet = null;
			WritableSheet xDensitySheet = null;
			WritableSheet yDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook( f );
				xHistogramSheet = XLSUtil.createNewPage( wb, "X histogram" );
				yHistogramSheet = XLSUtil.createNewPage( wb, "Y histogram" );
				xDensitySheet = XLSUtil.createNewPage( wb, "X density" );
				yDensitySheet = XLSUtil.createNewPage( wb, "Y density" );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
			XLSUtil.setCellString( xHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( xDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( xHistogramSheet, 0, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( xDensitySheet, 0, bin + 1, density1[ bin ] );
			}
			XLSUtil.setCellString( yHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( yDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( yHistogramSheet, 0, bin + 1, histogram2[ bin ] );
				XLSUtil.setCellNumber( yDensitySheet, 0, bin + 1, density2[ bin ] );
			}
			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
	}

	// export histograms and densities as xls file
	void exportSphericalFiles( final float[] histogram1, final float[] histogram2, final float[] histogram3, final float[] density1, final float[] density2, final float[] density3, final File f, final Sequence seq )
	{
		if ( f.exists() )
		{
			WritableWorkbook wb = null;
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angle1HistogramSheet = null;
			WritableSheet angle2HistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angle1DensitySheet = null;
			WritableSheet angle2DensitySheet = null;
			int radiusRows = 0, angle1Rows = 0, angle2Rows = 0, currentColumn = 0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite( f );
				if ( wb.getNumberOfSheets() != 6 )
				{
					MessageDialog.showDialog( "The excel output file should contain 6 worksheets: Radius histogram, Colatitude histogram, Azimuth angle histogram, Radius density, Colatitude density and Azimuth angle density." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				radiusHistogramSheet = wb.getSheet( 0 );
				angle1HistogramSheet = wb.getSheet( 1 );
				angle2HistogramSheet = wb.getSheet( 2 );
				radiusDensitySheet = wb.getSheet( 3 );
				angle1DensitySheet = wb.getSheet( 4 );
				angle2DensitySheet = wb.getSheet( 5 );
				currentColumn = radiusHistogramSheet.getColumns();
				if ( ( angle1HistogramSheet.getColumns() != currentColumn ) || ( angle2HistogramSheet.getColumns() != currentColumn ) || ( radiusDensitySheet.getColumns() != currentColumn ) || ( angle2DensitySheet.getColumns() != currentColumn ) || ( angle2DensitySheet.getColumns() != currentColumn ) )
				{
					MessageDialog.showDialog( "The excel output file should have the same number of columns in each worksheet." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				radiusRows = radiusHistogramSheet.getRows();
				if ( radiusDensitySheet.getRows() != radiusRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Radius histogram and Radius density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				angle1Rows = angle1HistogramSheet.getRows();
				if ( angle1DensitySheet.getRows() != angle1Rows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Colatitude histogram and Colatitude density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				angle2Rows = angle2HistogramSheet.getRows();
				if ( angle2DensitySheet.getRows() != angle2Rows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Azimuth angle histogram and Azimuth angle density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}

			if ( histogram1.length != ( radiusRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Radius histogram and Radius density do not have the same number of rows than the number of bins for radius." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( radiusHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( radiusDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( radiusHistogramSheet, currentColumn, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( radiusDensitySheet, currentColumn, bin + 1, density1[ bin ] );
			}
			if ( histogram2.length != ( angle1Rows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Colatitude histogram and Colatitude density do not have the same number of rows than the number of bins for colatitude." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( angle1HistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( angle1DensitySheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellNumber( angle1HistogramSheet, currentColumn, 1, histogram2[ 0 ] );
			XLSUtil.setCellNumber( angle1DensitySheet, currentColumn, 1, density2[ 0 ] );
			for ( int bin = 1; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( angle1HistogramSheet, currentColumn, bin + 1, histogram2[ histogram2.length - bin ] );
				XLSUtil.setCellNumber( angle1DensitySheet, currentColumn, bin + 1, density2[ histogram2.length - bin ] );
			}
			if ( histogram3.length != ( angle2Rows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Azimuth angle histogram and Azimuth angle density do not have the same number of rows than the number of bins for azimuth angle." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( angle2HistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( angle2DensitySheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellNumber( angle2HistogramSheet, currentColumn, 1, histogram3[ 0 ] );
			XLSUtil.setCellNumber( angle2DensitySheet, currentColumn, 1, density3[ 0 ] );
			for ( int bin = 1; bin < histogram3.length; bin++ )
			{
				XLSUtil.setCellNumber( angle2HistogramSheet, currentColumn, bin + 1, histogram2[ histogram3.length - bin ] );
				XLSUtil.setCellNumber( angle2DensitySheet, currentColumn, bin + 1, density2[ histogram3.length - bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
		else
		{
			WritableSheet radiusHistogramSheet = null;
			WritableSheet angle1HistogramSheet = null;
			WritableSheet angle2HistogramSheet = null;
			WritableSheet radiusDensitySheet = null;
			WritableSheet angle1DensitySheet = null;
			WritableSheet angle2DensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook( f );
				radiusHistogramSheet = XLSUtil.createNewPage( wb, "Radius histogram" );
				angle1HistogramSheet = XLSUtil.createNewPage( wb, "Azimuth angle histogram" );
				angle2HistogramSheet = XLSUtil.createNewPage( wb, "Colatitude histogram" );
				radiusDensitySheet = XLSUtil.createNewPage( wb, "Radius density" );
				angle1DensitySheet = XLSUtil.createNewPage( wb, "Azimuth angle density" );
				angle2DensitySheet = XLSUtil.createNewPage( wb, "Colatitude density" );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
			XLSUtil.setCellString( radiusHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( radiusDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram1.length; bin++ )
			{
				XLSUtil.setCellNumber( radiusHistogramSheet, 0, bin + 1, histogram1[ bin ] );
				XLSUtil.setCellNumber( radiusDensitySheet, 0, bin + 1, density1[ bin ] );
			}
			XLSUtil.setCellString( angle1HistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( angle1DensitySheet, 0, 0, seq.getName() );
			XLSUtil.setCellNumber( angle1HistogramSheet, 0, 1, histogram2[ 0 ] );
			XLSUtil.setCellNumber( angle1DensitySheet, 0, 1, density2[ 0 ] );
			for ( int bin = 1; bin < histogram2.length; bin++ )
			{
				XLSUtil.setCellNumber( angle1HistogramSheet, 0, bin + 1, histogram2[ histogram2.length - bin ] );
				XLSUtil.setCellNumber( angle1DensitySheet, 0, bin + 1, density2[ histogram2.length - bin ] );
			}
			XLSUtil.setCellString( angle2HistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( angle2DensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram3.length; bin++ )
			{
				XLSUtil.setCellNumber( angle2HistogramSheet, 0, bin + 1, histogram3[ bin ] );
				XLSUtil.setCellNumber( angle2DensitySheet, 0, bin + 1, density3[ bin ] );
			}
			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
	}

	// export histogram or density as xls file
	void exportDistanceFiles( final float[] histogram, final float[] density, final File f, final Sequence seq )
	{
		if ( f.exists() )
		{
			WritableWorkbook wb = null;
			WritableSheet distanceHistogramSheet = null;
			WritableSheet distanceDensitySheet = null;
			int distanceRows = 0, currentColumn = 0;

			try
			{
				wb = XLSUtil.loadWorkbookForWrite( f );
				if ( wb.getNumberOfSheets() != 2 )
				{
					MessageDialog.showDialog( "The excel output file should contain 2 worksheets: Distance to cell border histogram and Distance to cell border density." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				distanceHistogramSheet = wb.getSheet( 0 );
				distanceDensitySheet = wb.getSheet( 1 );
				currentColumn = distanceHistogramSheet.getColumns();
				if ( distanceDensitySheet.getColumns() != currentColumn )
				{
					MessageDialog.showDialog( "The excel output file should have the same number of columns in each worksheet." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
				distanceRows = distanceHistogramSheet.getRows();
				if ( distanceDensitySheet.getRows() != distanceRows )
				{
					MessageDialog.showDialog( "In the excel file, the worksheets Distance to cell border histogram and Distance to cell border density should have the same number of rows." );
					try
					{
						XLSUtil.saveAndClose( wb );
					}
					catch ( final Exception e )
					{
						throw new IcyHandledException( e.getMessage() );
					}
					return;
				}
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}

			if ( histogram.length != ( distanceRows - 1 ) )
			{
				MessageDialog.showDialog( "In the excel file, the worksheets Distance to cell border histogram and Distance to cell border density do not have the same number of rows than the number of bins for distance to cell border." );
				try
				{
					XLSUtil.saveAndClose( wb );
				}
				catch ( final Exception e )
				{
					throw new IcyHandledException( e.getMessage() );
				}
				return;
			}
			XLSUtil.setCellString( distanceHistogramSheet, currentColumn, 0, seq.getName() );
			XLSUtil.setCellString( distanceDensitySheet, currentColumn, 0, seq.getName() );
			for ( int bin = 0; bin < histogram.length; bin++ )
			{
				XLSUtil.setCellNumber( distanceHistogramSheet, currentColumn, bin + 1, histogram[ bin ] );
				XLSUtil.setCellNumber( distanceDensitySheet, currentColumn, bin + 1, density[ bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
		else
		{
			WritableSheet distanceHistogramSheet = null;
			WritableSheet distanceDensitySheet = null;
			WritableWorkbook wb = null;

			try
			{
				wb = XLSUtil.createWorkbook( f );
				distanceHistogramSheet = XLSUtil.createNewPage( wb, "Distance histogram" );
				distanceDensitySheet = XLSUtil.createNewPage( wb, "Distance density" );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
			XLSUtil.setCellString( distanceHistogramSheet, 0, 0, seq.getName() );
			XLSUtil.setCellString( distanceDensitySheet, 0, 0, seq.getName() );
			for ( int bin = 0; bin < histogram.length; bin++ )
			{
				XLSUtil.setCellNumber( distanceHistogramSheet, 0, bin + 1, histogram[ bin ] );
				XLSUtil.setCellNumber( distanceDensitySheet, 0, bin + 1, density[ bin ] );
			}

			try
			{
				XLSUtil.saveAndClose( wb );
			}
			catch ( final Exception e )
			{
				throw new IcyHandledException( e.getMessage() );
			}
		}
	}

	// compute CEMD
	float[][] CEMD( final float[][] density, final int nbExperiments, final int nbBins )
	{
		final float[][] output = new float[ nbExperiments ][ nbExperiments ];
		final float[][][] intermediateDensities = new float[ nbExperiments ][ nbBins ][ nbBins ];
		for ( int xp = 0; xp < nbExperiments; xp++ )
		{
			for ( int k = 0; k < nbBins; k++ )
			{
				intermediateDensities[ xp ][ k ][ 0 ] = density[ xp ][ k ];
				for ( int bin = 1; bin < nbBins; bin++ )
				{
					int currentBin = bin + k;
					if ( currentBin >= nbBins )
					{
						currentBin -= nbBins;
					}
					intermediateDensities[ xp ][ k ][ bin ] = intermediateDensities[ xp ][ k ][ bin - 1 ] + density[ xp ][ currentBin ];
				}
			}
		}

		for ( int i = 0; i < nbExperiments; i++ )
		{
			for ( int u = i + 1; u < nbExperiments; u++ )
			{
				float minOutput = 0;
				for ( int bin = 0; bin < nbBins; bin++ )
				{
					minOutput += Math.abs( intermediateDensities[ i ][ 0 ][ bin ] - intermediateDensities[ u ][ 0 ][ bin ] );
				}
				for ( int k = 1; k < nbBins; k++ )
				{
					float currentOutput = 0;
					for ( int bin = 0; bin < nbBins; bin++ )
					{
						currentOutput += Math.abs( intermediateDensities[ i ][ k ][ bin ] - intermediateDensities[ u ][ k ][ bin ] );
					}
					if ( currentOutput < minOutput )
					{
						minOutput = currentOutput;
					}
				}
				output[ i ][ u ] = minOutput;
				output[ u ][ i ] = output[ i ][ u ];
			}
		}
		return output;
	}

	// compute EMD
	float[][] EMD( final float[][] density, final int nbExperiments, final int nbBins )
	{
		final float[][] output = new float[ nbExperiments ][ nbExperiments ], cumulatedDensities = density;
		for ( int bin = 1; bin < nbBins; bin++ )
		{
			for ( int xp = 0; xp < nbExperiments; xp++ )
			{
				cumulatedDensities[ xp ][ bin ] += cumulatedDensities[ xp ][ bin - 1 ];
			}
		}
		for ( int i = 0; i < nbExperiments; i++ )
		{
			for ( int u = i + 1; u < nbExperiments; u++ )
			{
				for ( int bin = 0; bin < nbBins; bin++ )
				{
					output[ i ][ u ] += Math.abs( cumulatedDensities[ i ][ bin ] - cumulatedDensities[ u ][ bin ] );
				}
				output[ u ][ i ] = output[ i ][ u ];
			}
		}
		return output;
	}

	@Override
	protected void initialize()
	{
		// number of bins
		// Cartesian coordinate system
		final int nbBinsForX = 100;
		final int nbBinsForY = 100;
		final int nbBinsForCartesianDepth = 100;
		// cylindrical coordinate system
		final int nbBinsForCylindricalRadius = 100;
		final int nbBinsForAngle = 180;
		final int nbBinsForCylindricalDepth= 100;
		// spherical coordinate system
		final int nbBinsForSphericalRadius = 100;
		final int nbBinsForFirstAngle = 180;
		final int nbBinsForSecondAngle = 180;
		// distance distribution (distance to cell border for example)
		final int nbBinsForDistance = 100;

		// coordinate system center
		// Cartesian coordinate system
		final int ref1Xcartesian = -1;
		final int ref1Ycartesian = -1;
		// cylindrical coordinate system
		final int ref1Xcylindrical = -1;
		final int ref1Ycylindrical = -1;
		// spherical coordinate system
		final int ref1Xspherical = -1;
		final int ref1Yspherical= -1;
		final int ref1Zspherical = -1;

		// second point to define coordinate system reference direction
		// Cartesian coordinate system
		final int ref2Xcartesian = -1;
		final int ref2Ycartesian = -1;
		// cylindrical coordinate system
		final int ref2Xcylindrical = -1;
		final int ref2Ycylindrical = -1;
		// spherical coordinate system
		final int ref2Xspherical = -1;
		final int ref2Yspherical = -1;

		// coordinate system
		final String[] coordinateSystemPossibilities = {"Cylindrical","Spherical","Cartesian","Distance to cell border"};
		final String coordinateSystem = coordinateSystemPossibilities[0];

		// input data
		input;

		// input cell mask for cell normalization (and forbidden region if defined)
		final image cellMaskDensities = new image;
		final image forbiddenRegionDensities = new image;

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
		final int arraySize=inputImage.getSizeX()*inputImage.getSizeY(),
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
		final float[] intensitySumOverTime = new float[arraySize];
		for(int t=0;t<nbFrames;t++){
			final float[][] inputArray = Array2DUtil.arrayToFloatArray(inputImage.getDataXYZ(t,channelOfInterestDensities.getValue()), inputImage.isSignedDataType());
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


			final float maxRadius=0;
			// number of trajectories
			final int nbTrajectories;
			component1 = new float[nbTrajectories];
			component2 = new float[nbTrajectories];
			component3 = new float[nbTrajectories];
			weight = new float[nbTrajectories];
			distance1 = new float[nbTrajectories];
			distance2 = new float[nbTrajectories];
			distance2 = new float[nbTrajectories];

			// compute normalizing distances for each spatial point so it's done once for all
			final float[][] cellSegment = new float[depth][arraySize];
			final int[][] insideXCoord = new int[depth][arraySize];
			final int[][] insideYCoord = new int[depth][arraySize];
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

			final int cpt=0;
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

				final float minX=1000000,minY=100000,maxX=0,maxY=0;
				final int nbTrajectories;
				component1 = new float[nbTrajectories];
				component2 = new float[nbTrajectories];
				component3 = new float[nbTrajectories];
				weight = new float[nbTrajectories];
				distance1 = new float[nbTrajectories];
				distance2 = new float[nbTrajectories];

				// compute normalizing distance
				// compute normalizing distances for each spatial point so it's done once for all
				final float[][] cellSegment = new float[depth][arraySize];
				final int[][] insideXCoord = new int[depth][arraySize];
				final int[][] insideYCoord = new int[depth][arraySize];
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

				final int cpt=0;
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
				final float bandwidthForGaussianDistributionForX = compute_rule_of_thumb_bandwidth(component1);

				// compute x histograms
				final float[] interpolatedXHistogram = new float[nbBinsForX.getValue()],
						xHistogram = computeSymmetricalHistogram(component1,weight,distance1,nbBinsForX.getValue(),minX,maxX,interpolatedXHistogram);

				// compute Gaussian kernel
				final float[] GaussianKernelForX = computeSymmetricalGaussianKernel(nbBinsForX.getValue(),bandwidthForGaussianDistributionForX,maxX);

				// compute x density
				final float[] xDensity = computeDensity(nbBinsForX.getValue(),xHistogram,GaussianKernelForX);

				// y
				// compute density for y distribution for each experiment
				final float bandwidthForGaussianDistributionForY = compute_rule_of_thumb_bandwidth(component2);

				// compute x histograms
				final float[] interpolatedYHistogram = new float[nbBinsForY.getValue()],
						yHistogram = computeSymmetricalHistogram(component2,weight,distance1,nbBinsForY.getValue(),minY,maxY,interpolatedYHistogram);

				// compute Gaussian kernel
				final float[] GaussianKernelForY = computeSymmetricalGaussianKernel(nbBinsForY.getValue(),bandwidthForGaussianDistributionForY,maxY);

				// compute x density
				final float[] yDensity = computeDensity(nbBinsForY.getValue(),yHistogram,GaussianKernelForY);

				// depth
				// compute density for depth distribution for each experiment
				if(depth>1){
					// compute density for depth distribution for each experiment
					final float bandwidthForGaussianDistributionForDepth = compute_rule_of_thumb_bandwidth(component3);

					// compute depth histograms
					final float[] interpolatedDepthHistogram = new float[nbBinsForCartesianDepth.getValue()],
							depthHistogram = computeDepthHistogram(component3,weight,distance2,nbBinsForCartesianDepth.getValue(),depth-1,interpolatedDepthHistogram);

					// compute Gaussian kernel for each experiment
					final float[] GaussianKernelForDepth = computeGaussianKernel(nbBinsForCartesianDepth.getValue(),bandwidthForGaussianDistributionForDepth,depth-1);

					// compute depth density
					final float[] depthDensity = computeDensity(nbBinsForCartesianDepth.getValue(),depthHistogram,GaussianKernelForDepth);

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

					final float maxRadius=0;
					final int nbTrajectories;
					component1 = new float[nbTrajectories];
					component2 = new float[nbTrajectories];
					component3 = new float[nbTrajectories];
					weight = new float[nbTrajectories];
					distance1 = new float[nbTrajectories];
					distance2 = new float[nbTrajectories];

					// compute normalizing distance
					// compute normalizing distances for each spatial point so it's done once for all
					final float[][] cellSegment = new float[depth][arraySize];
					final int[][] insideXCoord = new int[depth][arraySize];
					final int[][] insideYCoord = new int[depth][arraySize];
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

					final int cpt=0;
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
					final float bandwidthForGaussianDistributionForRadius = compute_rule_of_thumb_bandwidth(component1);

					// compute radius histograms
					// To modify distanceToZborder -> distanceToPlaneBorder
					final float[] radiusHistogram = computeRadiusHistogram(component1,weight,distance1,nbBinsForSphericalRadius.getValue());

					// compute Gaussian kernel for each experiment
					final float[] GaussianKernelForRadius = computeGaussianKernel(nbBinsForSphericalRadius.getValue(),bandwidthForGaussianDistributionForRadius,maxRadius);

					// compute radius density
					final float[] radiusDensity = computeDensity(nbBinsForSphericalRadius.getValue(),radiusHistogram,GaussianKernelForRadius);


					// azimuthal angle
					// compute density for azimuthal angle distribution for each experiment
					final float concentrationForVonMisesDistributionOfAzimuthalAngle = compute_rule_of_thumb_concentration(component2);

					// compute azimuthal angle histograms
					final float[] azimuthalAngleHistogram = computeCircularHistogram(component2,weight,distance1,nbBinsForFirstAngle.getValue());

					// compute von Mises kernel for each experiment
					final float[] vonMisesKernelOfAzimuthalAngle = computeVonMisesKernel(nbBinsForFirstAngle.getValue(),concentrationForVonMisesDistributionOfAzimuthalAngle);

					// compute azimuthal angle density
					final float[] azimuthalAngleDensity = computeCircularDensity(nbBinsForFirstAngle.getValue(),azimuthalAngleHistogram,vonMisesKernelOfAzimuthalAngle);


					// polar angle
					// compute density for polar angle distribution for each experiment
					final float concentrationForVonMisesDistributionOfPolarAngle = compute_rule_of_thumb_concentration(component3);

					// compute polar angle histograms
					final float[] polarAngleHistogram = computeCircularHistogram(component3,weight,distance1,nbBinsForSecondAngle.getValue());

					// compute von Mises kernel for each experiment
					final float[] vonMisesKernelOfPolarAngle = computeVonMisesKernel(nbBinsForSecondAngle.getValue(),concentrationForVonMisesDistributionOfPolarAngle);

					// compute polar angle density
					final float[] polarAngleDensity = computeCircularDensity(nbBinsForSecondAngle.getValue(),polarAngleHistogram,vonMisesKernelOfPolarAngle);

					File f = exportSphericalExcelFile.getValue(true);
					if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
					exportSphericalFiles(radiusHistogram,azimuthalAngleHistogram,polarAngleHistogram,radiusDensity,azimuthalAngleDensity,polarAngleDensity,f,inputImage);
				}

				else{
					if(coordinateSystem.getValue()=="Distance to cell border"){

						float maxDistance=0;
						int nbTrajectories=0;
						for(int t=0;t<nbFrames;t++){
							final float[][] inputArray = Array2DUtil.arrayToFloatArray(inputImage.getDataXYZ(t,channelOfInterestDensities.getValue()), inputImage.isSignedDataType());
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
						final float[][] cellSegment = new float[depth][arraySize];
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
						final float bandwidthForGaussianDistributionForDistance = compute_rule_of_thumb_bandwidth(component1);

						// compute distance histograms
						final float[] distanceHistogram = computeRadiusHistogram(component1,weight,distance1,nbBinsForDistance.getValue());

						// compute Gaussian kernel for each experiment
						final float[] GaussianKernelForDistance = computeGaussianKernel(nbBinsForDistance.getValue(),bandwidthForGaussianDistributionForDistance,maxDistance);

						// compute distance density
						final float[] distanceDensity = computeDensity(nbBinsForDistance.getValue(),distanceHistogram,GaussianKernelForDistance);

						File f = exportDistanceExcelFile.getValue(true);
						if (!FileUtil.getFileExtension(f.getPath(), false).equalsIgnoreCase("xls")) f = new File(f.getPath() + ".xls");
						exportDistanceFiles(distanceHistogram,distanceDensity,f,inputImage);
					}
				}
			}

		}
	}

}