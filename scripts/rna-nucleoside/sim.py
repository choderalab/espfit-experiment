import click
from espfit.app.sampler import SetupSampler


def run(kwargs):
    # Parameters
    configfile = kwargs['configfile']
    override_sampler_kwargs = {
        'override_with_espaloma': kwargs['override_with_espaloma'],
        'small_molecule_forcefield': kwargs['small_molecule_forcefield'],
    }

    # Setup
    sampler = SetupSampler._from_toml(configfile, **override_sampler_kwargs)[0]
    sampler.export_xml(exportSystem=True, exportState=False, exportIntegrator=True)
    
    # Debug
    print(f'small molecule forcefield: {sampler.small_molecule_forcefield}')
    print(f'override_with_espaloma: {sampler.override_with_espaloma}')
    print(f'forcefield list: {sampler.forcefield_files}')
    print(f'water model: {sampler.water_model}, water class: {sampler.water_class}')
    print(f'maxIterations: {sampler.maxIterations}')
    print(f'nsteps: {sampler.nsteps}')
    
    # Minimize and run
    sampler.minimize()
    sampler.run()

    # Export
    sampler.export_xml(exportSystem=False, exportState=True, exportIntegrator=False)


@click.command()
@click.option('--configfile', default='config.toml', help='Configuration toml file', type=str)
@click.option('--small_molecule_forcefield', default='openff-1.2.1.offxml', help='Forcefield for small molecule', type=str)
@click.option('--override_with_espaloma', default=True, help='Override the solute system with espaloma', type=bool)
def cli(**kwargs):
    run(kwargs)


if __name__ == "__main__":
    cli()