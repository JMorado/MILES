Testing plugins (provisional)
#############################

Fixtures and helpers for a plugin's tests: ``pytest_plugins = ["miles.testing"]`` in the plugin's ``conftest.py``, then ``from miles.testing import ...``. :ref:`Testing a plugin <dg-plugin-testing>` shows them at work. Provisional: they follow the window of MILES and may change between releases (see :ref:`dg-stability`).

.. automodule:: miles.testing
   :no-members:

Fixtures
********

.. py:function:: miles.testing.miles_app()
   :no-index:

   A new MILES window with an empty tab, as an :class:`~miles.plugin.App`; it has settings of its own and is closed after the test without asking to save.

.. py:function:: miles.testing.qt_application()
   :no-index:

   The Qt application, made once for all the tests.

The window
**********

.. autofunction:: miles.testing.load_plugin

.. autofunction:: miles.testing.document_from_smiles

.. autofunction:: miles.testing.document_from_text

.. autofunction:: miles.testing.wait_until

The plugin's life
*****************

.. autofunction:: miles.testing.switch_off

.. autofunction:: miles.testing.switch_on

.. autofunction:: miles.testing.reload_plugin

.. autofunction:: miles.testing.plugin_setting

What the user sees
******************

.. autofunction:: miles.testing.captured_messages

.. autoclass:: miles.testing.Message

.. autofunction:: miles.testing.drawings
