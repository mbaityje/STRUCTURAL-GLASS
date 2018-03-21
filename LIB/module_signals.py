import signal
class GracefulKiller:
	'''
	Class to catch the sigterm,sigint signals and exit after having performed the due operations.
	'''
	kill_now = False
	def __init__(self):
		''' Custom signal handling '''
		signal.signal(signal.SIGINT, self.exit_gracefully)
		signal.signal(signal.SIGTERM, self.exit_gracefully)
		return

	def resetSignals(self):
		''' Reset default signal handling '''
		signal.signal(signal.SIGINT, signal.SIG_DFL)
		signal.signal(signal.SIGTERM, signal.SIG_DFL)
		return

	def exit_gracefully(self,signum, frame):
		self.kill_now = True
		return
