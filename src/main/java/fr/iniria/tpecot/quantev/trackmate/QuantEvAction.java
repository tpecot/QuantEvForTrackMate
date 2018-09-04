package fr.iniria.tpecot.quantev.trackmate;

import javax.swing.ImageIcon;

import org.scijava.plugin.Plugin;

import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.TrackModel;
import fiji.plugin.trackmate.action.AbstractTMAction;
import fiji.plugin.trackmate.action.CaptureOverlayAction;
import fiji.plugin.trackmate.action.TrackMateAction;
import fiji.plugin.trackmate.action.TrackMateActionFactory;
import fiji.plugin.trackmate.gui.TrackMateGUIController;
import fiji.plugin.trackmate.visualization.trackscheme.TrackSchemeFrame;
import ij.ImagePlus;

public class QuantEvAction extends AbstractTMAction implements TrackMateAction
{

	public static final ImageIcon ICON = new ImageIcon( TrackSchemeFrame.class.getResource( "QuantEv_pluginIcon.png" ) );

	public static final String NAME = "QuantEv";

	public static final String KEY = "QUANTEV";

	public static final String INFO_TEXT = "<html>" +
			"Runs the QuantEV analysis."
			+ "<p>"
			+ "See Pécot <i>et al.</i>: A quantitative approach for analyzing the spatio-temporal distribution of "
			+ "3D intracellular events in fluorescence microscopy. "
			+ "<code>https://elifesciences.org/articles/32311</code>" +
			"</html>";


	@Override
	public void execute( final TrackMate trackmate )
	{
		final TrackModel trackModel = trackmate.getModel().getTrackModel();
		final int nTracks = trackModel.nTracks( true );

		final ImagePlus imp = trackmate.getSettings().imp;
		if (imp == null)
		{
			logger.error( "Cannot run QuantEv: The target image is not set." );
			return;
		}

		// TODO
		System.out.println( "TODO: Run QuantEv on image " +imp.getShortTitle() + " with " + nTracks + " tracks." ); // TODO
	}

	@Plugin( type = TrackMateActionFactory.class )
	public static class Factory implements TrackMateActionFactory
	{

		@Override
		public String getInfoText()
		{
			return INFO_TEXT;
		}

		@Override
		public String getKey()
		{
			return KEY;
		}

		@Override
		public TrackMateAction create( final TrackMateGUIController controller )
		{
			return new CaptureOverlayAction( controller.getGUI() );
		}

		@Override
		public ImageIcon getIcon()
		{
			return ICON;
		}

		@Override
		public String getName()
		{
			return NAME;
		}

	}

}
