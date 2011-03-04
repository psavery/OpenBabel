/**********************************************************************
  SSHManager - Manages a collection of SSHConnections

  Copyright (C) 2010 by David C. Lonie

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#ifndef SSHMANAGER_H
#define SSHMANAGER_H

#include <globalsearch/sshconnection.h>

#include <QtCore/QObject>
#include <QtCore/QMutex>
#include <QtCore/QHash>

namespace GlobalSearch {
  class OptBase;

  /**
   * @class SSHManager sshmanager.h <globalsearch/sshmanager.h>
   *
   * @brief A class to manage multiple SSHConnection objects.
   *
   * @author David C. Lonie
   */
  class SSHManager : public QObject
  {
    Q_OBJECT

  public:
    /**
     * Constructor.
     *
     * @param connections The maximum number of simultaneous connections.
     * @param parent The OptBase parent
     */
    explicit SSHManager(unsigned int connections = 5, OptBase *parent = 0);

    /**
     * Destructor.
     */
    virtual ~SSHManager();

    /**
     * Create connections to the specifed host. If the connections
     * cannot be made, an SSHConnection::SSHConnectionException will
     * be thrown.
     */
    void makeConnections(const QString &host,
                         const QString &user = "",
                         const QString &pass = "",
                         unsigned int port = 22);

    /**
     * @return Whether the connection has been made successfully.
     */
    bool isValid() {return m_isValid;};

    /// Get the currently set user name
    QString getUser() {return m_user;};

    /// Get the currently set hostname
    QString getHost() {return m_host;};

    /// Get the currently set port
    int getPort() {return m_port;};

  public slots:
    /**
     * Returns a free connection from the pool and locks it.
     * @sa unlockConnection
     */
    SSHConnection *getFreeConnection();

    /**
     * Call this when finished with a connection so other threads can
     * use it.
     */
    void unlockConnection(SSHConnection* ssh);

    /**
     * Retreive the public key from the server. This is set when a
     * connection fails with SSH_UNKNOWN_HOST_ERROR.
     *
     * @sa SSH_UNKNOWN_HOST_ERROR
     * @sa validateServerKey
     */
    QString getServerKeyHash();

    /**
     * Add currently set key to the known host cache.
     *
     * @sa SSH_UNKNOWN_HOST_ERROR
     * @sa getServerKey;
     */
    bool validateServerKey();

    /**
     * Set the server key. This is used internally.
     */
    void setServerKey(const QString &hexa);


  protected:
    /// List of all SSHConnection objects managed by this instance
    QList<SSHConnection*> m_conns;

    /// Internally used mutex
    QMutex m_lock;

    /// Hostname
    QString m_host;
    /// Username
    QString m_user;
    /// Password
    QString m_pass;
    /// Key
    QString m_hexa;
    /// Port
    unsigned int m_port;
    /// Number of connections
    unsigned int m_connections;
    /// Monitor whether the connections are valid
    bool m_isValid;
  };

} // end namespace GlobalSearch

#endif