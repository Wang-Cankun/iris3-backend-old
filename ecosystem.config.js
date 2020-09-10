module.exports = {
  apps: [
    {
      name: 'IRIS3-backend',
      script: 'dist/main.js',
      // pm2 logs 0 --lines 1000
      // Options reference: https://pm2.io/doc/en/runtime/reference/ecosystem-file/
      instances: 1,
      autorestart: true,
      watch: ['.env'],
      ignore_watch: ['node_modules'],
      log_date_format: 'YYYY-MM=DD HH:mm:ss',
      max_memory_restart: '16G',
      env: {
        NODE_ENV: 'production'
      }
    }
  ]
}
