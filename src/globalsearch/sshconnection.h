/**********************************************************************
  SSHConnection - Connection to an ssh server for execution, sftp, etc.

  Copyright (C) 2010 by David C. Lonie

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#ifndef SSHCONNECTION_H
#define SSHCONNECTION_H

#define LIBSSH_STATIC
extern "C" {
#include <libssh/libssh/libssh.h>
#include <libssh/libssh/sftp.h>
}

#include <globalsearch/optbase.h>

#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtCore/QMutex>

#define SSH_BUFFER_SIZE 20480

namespace GlobalSearch {

  /**
   * @class SSHConnection sshconnection.h <globalsearch/sshconnection.h>
   *
   * @brief A class to handle command execution and sftp transactions
   * on an ssh server.
   *
   * This class should not be created directly. Instead, see SSHManager.
   *
   * @author David C. Lonie
   */
  class SSHConnection : public QObject
  {
    Q_OBJECT

  public:
    /// Exceptions
    enum SSHConnectionException {
      /// An error connecting to the host has occurred
      SSH_CONNECTION_ERROR = 1,
      /// The host is unknown or has changed its key.
      /// The key can be retrieved through SSHManager::getServerKey()
      /// and accepted via SSHManager::validateServerKey()
      SSH_UNKNOWN_HOST_ERROR,
      /// A bad password was given and public key auth is not set up
      SSH_BAD_PASSWORD_ERROR,
      /// An unknown error has occurred
      SSH_UNKNOWN_ERROR
    };

    /**
     * Constructor.
     *
     * @param parent The OptBase parent
     */
    explicit SSHConnection(SSHManager *parent = 0);

    /**
     * Destructor.
     */
    virtual ~SSHConnection();

    /// @return The currently set username
    QString getUser() {return m_user;};
    /// @return The currently set hostname
    QString getHost() {return m_host;};
    /// @return The currently set port
    int getPort() {return m_port;};

  public slots:
    /**
     * Set the login details for this connections
     *
     * @param host Hostname
     * @param user Username
     * @param pass Password
     * @param port Port
     */
    void setLoginDetails(const QString &host,
                         const QString &user = "",
                         const QString &pass = "",
                         int port = 22);

    /// Flag whether this connection is in use or not.
    void setUsed(bool b) {m_inUse = b;};

    /// @return True if this connection is in use.
    bool inUse() {return m_inUse;};

    /**
     * Execute an arbitrary command on the connected host.
     *
     * @param command Command to execute
     * @param stdout_str (return) standard output
     * @param stderr_str (return) standard output
     * @param exitcode (return) exit code.
     *
     * @return True on success.
     */
    bool execute(const QString &command,
                 QString &stdout_str,
                 QString &stderr_str,
                 int &exitcode);

    /**
     * Copy a file to the remote host
     *
     * @param localpath Local file to copy
     * @param remotepath Destination path on connected host
     *
     * @return True on success.
     */
    bool copyFileToServer(const QString & localpath,
                          const QString & remotepath);

    /**
     * Copy a file from the remote host
     *
     * @param remotepath Remote file to copy
     * @param localpath Destination path on local host
     *
     * @return True on success.
     */
    bool copyFileFromServer(const QString & remotepath,
                            const QString & localpath);

    /**
     * Obtain a QString object contain the contents of a remote text
     * file.
     *
     * @param filename Name of text file on remote host to read
     * @param contents (return) contents of file.
     *
     * @return True on success.
     */
    bool readRemoteFile(const QString &filename,
                        QString &contents);

    /**
     * Delete a file from the remote host
     *
     * @param filename Remote file to remove.
     *
     * @return True on success.
     */
    bool removeRemoteFile(const QString &filename);

    /**
     * Copy a directory to the remote host
     *
     * @param localpath Local directory to copy
     * @param remotepath Destination path on connected host
     *
     * @return True on success.
     */
    bool copyDirectoryToServer(const QString & localpath,
                               const QString & remotepath);

    /**
     * Copy a directory from the remote host
     *
     * @param remotepath Remote directory to copy
     * @param localpath Destination path on local host
     *
     * @return True on success.
     */
    bool copyDirectoryFromServer(const QString & remotepath,
                                 const QString & localpath);

    /**
     * List the contents of a remote directory
     *
     * @param remotepath Path to remote directory
     * @param contents (return) List of file/directory names
     *
     * @return True on success.
     */
    bool readRemoteDirectoryContents(const QString & remotepath,
                                     QStringList & contents);

    /**
     * Recursively delete a directory and its contents.
     *
     * @param remotepath Directory on remote host to delete
     * @param onlyDeleteContents If true, only clean the directory
     * (default: false)
     *
     * @return True on success.
     */
    bool removeRemoteDirectory(const QString & remotepath,
                               bool onlyDeleteContents = false);


    /// @return True if the session is valid
    bool isValid() {return m_isValid;};

    /// @return True if the session is connected
    bool isConnected();

    /**
     * Attempts to create a connection to the remote server. If \a
     * throwExceptions is true, one of the
     * GlobalSearch::SSHConnection::SSHConnectionException exceptions
     * will be thrown on errors. Otherwise, the function will return
     * false on errors.
     *
     * @return True on success.
     */
    bool connectSession(bool throwExceptions = false);

    /**
     * Disconnect and connect this session.
     *
     * @sa disconnectSession
     * @sa connectSession
     */
    bool reconnectSession(bool throwExceptions = false);

    /**
     * Terminate the current connection to the remote host.
     *
     * @return True on success.
     */
    bool disconnectSession();

    /**
     * Similar to reconnectSession, but calls isConnected first and
     * only reconects if needed.
     *
     * @return True if no errors.
     */
    bool reconnectIfNeeded() {if (!isConnected()) return reconnectSession(false);
      return true;};

    /**
     * Add the passed host's key to the local host's knownhosts cache.
     *
     * @param host Hostname of server
     * @param port Port of SSH server
     *
     * @return True on success.
     */
    static bool addKeyToKnownHosts(const QString &host, unsigned int port = 22);

  signals:
    /**
     * Emitted when a host is not recognized
     *
     * @param hexa Hex data representing the server.
     */
    void unknownHostKey(const QString &hexa);

  protected:
    // Disable doxygen parsing
    /// \cond
    bool _execute(const QString &command,
                  QString &stdout_err,
                  QString &stderr_err,
                  int &exitcode);
    bool _copyFileToServer(const QString & localpath,
                           const QString & remotepath);
    bool _copyFileFromServer(const QString & remotepath,
                             const QString & localpath);
    bool _readRemoteFile(const QString &filename,
                         QString &contents);
    bool _removeRemoteFile(const QString &filename);
    bool _copyDirectoryToServer(const QString & localpath,
                                const QString & remotepath);
    bool _copyDirectoryFromServer(const QString & remotepath,
                                  const QString & localpath);
    bool _readRemoteDirectoryContents(const QString & remotepath,
                                      QStringList & contents);
    bool _removeRemoteDirectory(const QString & remotepath,
                                bool onlyDeleteContents = false);

    ssh_session m_session;
    ssh_channel m_shell;
    sftp_session m_sftp;
    QString m_host;
    QString m_user;
    QString m_pass;
    int m_port;
    bool m_isValid;
    bool m_inUse;
    QMutex m_lock;

    // Resume doxygen parsing
    /// \endcond
  };

} // end namespace GlobalSearch

#endif