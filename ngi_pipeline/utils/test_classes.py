from ngi_pipeline.utils.classes import with_ngi_config

# This isn't being called by nosetests, I think due to the fact
# that it gets renamed as "with_config" despite the fact that I've
# attempted to address this inside the with_config definition
# with functools.update_wrapper. Wtf!
@with_ngi_config
def test_with_ngi_config_no_path(config=None, config_file_path=None):
    assert(config)
